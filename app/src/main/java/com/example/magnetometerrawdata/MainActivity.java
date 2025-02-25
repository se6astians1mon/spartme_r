package com.example.magnetometerrawdata;

import android.content.Context;
import android.content.Intent;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.os.Bundle;
import android.os.CountDownTimer;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import androidx.appcompat.app.AppCompatActivity;

import com.google.android.material.tabs.TabLayout;
import com.opencsv.CSVReader;

import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

public class MainActivity extends AppCompatActivity {
    private TextView zValueField;
    private TextView yValueField;
    private TextView xValueField;
    private TextView zValueAvgField;
    private TextView yValueAvgField;
    private TextView xValueAvgField;
    private EditText x;
    private EditText y;
    private TextView realx;
    private TextView realy;
    private TabLayout tabLayout;
    private Context context;
    private Button btn_calibrate;
    private SensorManager sensorManager;
    private SensorEventListener magnetSensorListener;
    private Sensor mag;
    private boolean record;
    private CountDownTimer timer;
    private int realX;
    private int realY;
    Queue<Float> samplesX;
    Queue<Float> samplesY;
    Queue<Float> samplesZ;
    int bufferLength;
    float[] calibrationValues;
    CustomView c;
    boolean log = false;
    Button btn_log;
    Queue<double[]> bufferDirectionVector;
    Queue<float[]> bufferX;
    private int bufferLengthDisplay;

    float a = 0, b = 0;
    float y1 = 0;
    private float x1 = 0, x2 = 0, y2 = 0;
    Button btn_setVal;
    EditText edt_cmdistance;
    TextView txt_cmdistance;
    TextView lbl_Position;
    private double length;
    private Queue<Float> bufferPosX;
    private Queue<Float> bufferPosY;


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        SensorManager sensorManager = (SensorManager) getSystemService(SENSOR_SERVICE);
        //switch to know when to record data and send it and when only to display new values
        record = false;
        samplesX = new LinkedList<Float>();
        samplesY = new LinkedList<Float>();
        samplesZ = new LinkedList<Float>();
        bufferDirectionVector = new LinkedList<double[]>();

        bufferX = new LinkedList<float[]>();
        bufferPosX = new LinkedList<>();
        bufferPosY = new LinkedList<>();
        bufferLengthDisplay = 400;
        bufferLength = 10;
        calibrationValues = new float[]{0.0f, 0.0f, 0.0f};
        //init queue
        for(int i = 0;i<bufferLength;i++){
            samplesX.add(0.0f);
            samplesY.add(0.0f);
            samplesZ.add(0.0f);
            bufferDirectionVector.add(new double[]{0.0f, 0.0f});
            bufferPosX.add(0.0f);
            bufferPosY.add(0.0f);
        }
        for(int i = 0;i<bufferLengthDisplay;i++){
            bufferX.add(new float[]{0.0f, 0.0f});
            bufferDirectionVector.add(new double[]{0.0f, 0.0f});

        }
        magnetSensorListener = new SensorEventListener() {
            @Override
            public void onSensorChanged(SensorEvent event) {
                if (event.sensor.getType() == Sensor.TYPE_MAGNETIC_FIELD) {
                    // get values for each axes X,Y,Z
                    zValueField.setText(Float.toString(event.values[2]-calibrationValues[2]));
                    xValueField.setText(Float.toString(event.values[0]-calibrationValues[0]));
                    yValueField.setText(Float.toString(event.values[1]-calibrationValues[1]));

                    //calc running sum
                    samplesX.remove();
                    samplesY.remove();
                    samplesZ.remove();
                    samplesX.add(event.values[0]);
                    samplesY.add(event.values[1]);
                    samplesZ.add(event.values[2]);
                    if(log)
                        Log.d("sensorvalues",Float.toString(avg(samplesX)-calibrationValues[0])+","+Float.toString(avg(samplesY)-calibrationValues[1]));
                    xValueAvgField.setText(Float.toString(avg(samplesX)-calibrationValues[0]));
                    yValueAvgField.setText(Float.toString(avg(samplesY)-calibrationValues[1]));
                    zValueAvgField.setText(Float.toString(avg(samplesZ)-calibrationValues[2]));

                    bufferX.remove();
                    bufferX.add(new float[]{avg(samplesX)-calibrationValues[0], avg(samplesY)-calibrationValues[1]});
                    //try to guess ellipse:





                    int[] start_end_ellipse = findEllipseThroughAutocorrelation(bufferX, true);
                    //Log.i("eclipsefinding", start_end_ellipse[0]+" ,"+start_end_ellipse[1]);
                    if(start_end_ellipse[1]>0){
                        double[] xX = new double[bufferX.size()];
                        double[] yY = new double[bufferX.size()];
                        int i = 0;
                        for (float[] point : bufferX) {
                            xX[i] = point[0];
                            yY[i] = point[1];
                            i++;
                        }
                        try {
                            double[] params = EllipseFitter.fitEllipse(xX, yY);
                            double[] fittedParams = EllipseFitter.cartToPol(params);
                            Log.i("directellipsefit", "Fitted Ellipse Parameters: " + fittedParams[0]+", "+fittedParams[1]);
                        } catch (IllegalArgumentException e) {
                            Log.e("Ellipse Error"," "+e.getMessage());
                        }

                        double[][] points = PCAEllipse.queueTo2DArray(bufferX, start_end_ellipse);
                        double[] directionvector =  rotateVector(PCAEllipse.findMajor(points), 0);
                        c.center = PCAEllipse.getCenter(points);
                        //find which side of the model the stuff is largest away (measure distance between points and max)
                        List<double[]> pointsList = new ArrayList<>();
                        for (float[] point : bufferX) {
                            pointsList.add(new double[]{point[0], point[1]});
                        }
                        int indexMaxSpace = findMaxSpacing(pointsList, 10);
                        double[] widestPoint = pointsList.get(indexMaxSpace);

                        c.pointSpace[0] = (float)widestPoint[0];
                        c.pointSpace[1] = (float)widestPoint[1];
                        //if left, make sure to flip vector
                        //translate widest point to center...
                        widestPoint[0] -= c.center[0];
                        widestPoint[1] -= c.center[1];
                        if(pointDirection(avgVector(bufferDirectionVector), widestPoint).equals("left")){
                            Log.i("direction", "to my left I need to switch!!");
                            c.vector[0] = -avgVector(bufferDirectionVector)[0];
                            c.vector[1] = -avgVector(bufferDirectionVector)[1];
                        } else {
                            c.vector = avgVector(bufferDirectionVector);
                            Log.i("direction", "to my right, everything is al right!");
                        }
                        float angle = vector2angle(c.vector);



                        bufferDirectionVector.add(new double[]{directionvector[0], directionvector[1]});
                        bufferDirectionVector.remove();

                        length = PCAEllipse.calculateMajorAxisLength(points);

                        c.length = (float)length;
                        realx.setText("Length [µT]: " +String.format("%.2f",length)+"");
                        double distancecm = Math.sqrt(a/(length-b));
                        double distancecmModel = 120.413*Math.pow(length, -0.352761);
                        txt_cmdistance.setText(distancecmModel+" cm");
                        Log.i("trig",""+angle+","+distancecmModel);
                        double x = distancecmModel * Math.cos(angle);
                        double y = distancecmModel * Math.sin(angle);
                        if(distancecmModel<100.0)
                            Log.i("pos",String.format("%.1f",avg(bufferPosX))+"|"+String.format("%.1f",avg(bufferPosY))+"|"+String.format("%.1f",distancecmModel)+"|"+String.format("%.1f",angle));
                        bufferPosY.remove();
                        bufferPosY.add(new Float(y));
                        bufferPosX.remove();
                        bufferPosX.add(new Float(x));
                        lbl_Position.setText(String.format("%.1f",avg(bufferPosX))+","+String.format("%.1f",avg(bufferPosY)));


                    }

                    c.data = new float[]{avg(samplesX)-calibrationValues[0], avg(samplesY)-calibrationValues[1], avg(samplesZ)-calibrationValues[2]};
                    c.bufferX = bufferX;
                    c.bufferLength = bufferLengthDisplay;
                    float a = avg(samplesX)-calibrationValues[0];
                    float b = avg(samplesY)-calibrationValues[1];
                    double c_ = Math.sqrt(a*a+b*b);

                    c.invalidate();
                }
            }

            @Override
            public void onAccuracyChanged(Sensor sensor, int i) {

            }
        };
        mag = sensorManager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD_UNCALIBRATED);

        //init all the gui refs.
        c = findViewById(R.id.customView);
        realx = (TextView) findViewById(R.id.txt_xyvectorlength);
        realy = (TextView) findViewById(R.id.txt_xyvectorlength);
        lbl_Position = findViewById(R.id.lbl_Position);

        zValueField = (TextView) findViewById(R.id.txt_z);
        xValueField = (TextView) findViewById(R.id.txt_x);
        yValueField = (TextView) findViewById(R.id.txt_y);

        zValueAvgField = (TextView) findViewById(R.id.txt_z_avg);
        xValueAvgField = (TextView) findViewById(R.id.txt_x_avg);
        yValueAvgField = (TextView) findViewById(R.id.txt_y_avg);
        edt_cmdistance = findViewById(R.id.editTextNumber);
        txt_cmdistance = findViewById(R.id.txt_cmdistance);
        btn_setVal = findViewById(R.id.btn_setDistance);
        btn_setVal.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(btn_setVal.getText().equals("Map Distance")){
                    y1 = Float.parseFloat(edt_cmdistance.getText().toString());
                    x1 = (float) length;
                    btn_setVal.setText("Set 2nd Value");

                } else if(btn_setVal.getText().equals("Set 2nd Value")){
                    btn_setVal.setText("Map Distance");
                    y2 = Float.parseFloat(edt_cmdistance.getText().toString());
                    x2 = (float) length;
                    //find values for f(x) = a/x²+b
                    b =  ((y2*x2*x2)-(y1*x1*x1))/(x2*x2-x1*x1);
                    a = (y1-b)*x1*x1;
                }
            }
        });


        context = this;
        tabLayout = findViewById(R.id.tabLayout);
        tabLayout.getTabAt(0).select();
        tabLayout.addOnTabSelectedListener(new TabLayout.OnTabSelectedListener() {
            @Override
            public void onTabSelected(TabLayout.Tab tab) {
                Log.i("info", "chouette " + tab.getText());
                if(tab.getText().equals("prototype")){
                    Log.i("info", "chouette equals acc");
                    Intent intent = new Intent(context, PrototypeActivity.class);
                    startActivity(intent);
                }
            }

            @Override
            public void onTabUnselected(TabLayout.Tab tab) {

            }

            @Override
            public void onTabReselected(TabLayout.Tab tab) {

            }
        });
        btn_log = findViewById(R.id.btn_log);
        btn_log.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                log = !log;
                btn_log.setText("logging: "+log );
            }
        });
        btn_calibrate = findViewById(R.id.btn_calibrate);
        btn_calibrate.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                calibrationValues[0]= avg(samplesX);
                calibrationValues[1]= avg(samplesY);
                calibrationValues[2]= avg(samplesZ);
                btn_calibrate.setText("Calibrated");
            }
        });
        sensorManager.registerListener(magnetSensorListener, sensorManager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD),100);

        //get mapping
        ArrayList<Float[]> mapping = importMapping();
    }
    //returns angle as degrees
    private float vector2angle(double[] vector) {
        double angleRadians = Math.atan2(vector[0], vector[1]);
        // Optionally, convert radians to degrees
        double angleDegrees = Math.toDegrees(angleRadians);
        angleDegrees = (angleDegrees + 360) % 360;
        return (float)angleDegrees; // Return degrees if needed
    }

    private double[] avgVector(Queue<double[]> bufferDirectionVector) {
        double[] avgVector = new double[]{0.0f, 0.0f};
        for (double[] vector :
                bufferDirectionVector) {
            avgVector[0] += vector[0];
            avgVector[1] += vector[1];
        }
        avgVector[0]/=bufferDirectionVector.size();
        avgVector[1]/=bufferDirectionVector.size();
        return avgVector;
    }

    // returns position of first two lags in data. if none found returns 0 0
    private int[] findEllipseThroughAutocorrelation(Queue<float[]> data,boolean first) {
        double[] dataX = extractOneDimension(data, 0);
        double[] autocorr = computeAutocorrelation(dataX);
        List<Integer> maxima = findLocalMaxima(autocorr);

        // Step 3: Use the first two maxima to extract the ellipse data
        if (maxima.size() >= 2) {
            if (first){
                int start = maxima.get(0);
                int end = maxima.get(1);
                return new int[]{start, end};
            } else {
                int start = maxima.get(maxima.size()-2);
                int end = maxima.get(maxima.size()-1);
                return new int[]{start, end};
            }


            // Now you can process ellipseData further (e.g., finding the major axis)
        }
        else
            return new int[]{0,0};
    }

    private double[] extractOneDimension(Queue<float[]> data, int index) {
        int size = data.size();
        double[] values = new double[size];  // 2 columns for x and y
        int i = 0;
        for (float[] point : data) {
            values[i] = point[index];  // x value
            i++;
        }

        return values;  // Returns 2D array of (x, y) pairs
    }

    public static double[] computeAutocorrelation(double[] data) {
        int N = data.length;
        double[] autocorr = new double[N];
        double mean = calculateMean(data);

        // Compute denominator (variance of the data)
        double denom = 0.0;
        for (double datum : data) {
            denom += (datum - mean) * (datum - mean);
        }

        // Compute autocorrelation for each lag
        for (int tau = 0; tau < N; tau++) {
            double num = 0.0;
            for (int t = 0; t < N - tau; t++) {
                num += (data[t] - mean) * (data[t + tau] - mean);
            }
            autocorr[tau] = num / denom;
        }

        return autocorr;
    }

    // Function to find local maxima in the autocorrelation values
    public static List<Integer> findLocalMaxima(double[] autocorr) {
        List<Integer> maxima = new ArrayList<>();
        for (int i = 1; i < autocorr.length - 1; i++) {
            if (autocorr[i] > autocorr[i - 1] && autocorr[i] > autocorr[i + 1]) {
                maxima.add(i); // Store the index of the local maximum
            }
        }
        return maxima;
    }

    // Helper function to calculate mean
    private static double calculateMean(double[] data) {
        double sum = 0.0;
        for (double datum : data) {
            sum += datum;
        }
        return sum / data.length;
    }
    ArrayList<Float[]> importMapping(){
        //loads a mapping from real coordinates x [0] y [1] to magnetic field strength values x [2] y [3]
        ArrayList<Float[]> result = new ArrayList<>();
        try {
            Log.i("Info", "test");
            CSVReader reader = new CSVReader(new InputStreamReader(getResources().openRawResource(R.raw.output)));
            String[] nextLine;
            boolean firstLine=true;
            while ((nextLine = reader.readNext()) != null) {
                // nextLine[] is an array of values from the line
                if(firstLine){
                    firstLine=false;
                    continue;
                }
                String xReal = nextLine[0].toString();
                String yReal = nextLine[1].toString();
                String xValue = nextLine[2].toString();
                String yValue = nextLine[3].toString();
                //Log.i("info", "values: "+xValue + ", "+yValue+", "+xReal+", "+yReal);
                result.add(new Float[] {Float.parseFloat(xReal),Float.parseFloat(yReal), Float.parseFloat(xValue), Float.parseFloat(yValue)});
            }
        } catch (IOException e) {

        }
        Log.i("Info", "imported "+result.size()+ " rows");
        return result;
    }
    private float avg(Queue<Float> samples) {
        float sum = 0.0f;
        for(float number : samples){
            sum+=number;
        }
        return sum/samples.size();
    }
    private float[] calcPos(ArrayList<Float[]> posList, float xt, float yt){
        //lookup: calc distance to each point in list. If distance < previous distance, overwrite previous distance with distance and move on, else move on. At the end return the real coords
        double previous_distance = 10000;
        float[] result = {0f,0f};
        float realX = 0f;
        float realY = 0f;
        for (Float[] f : posList){
            double distance = Math.sqrt((f[2]-xt)*(f[2]-xt)+(f[3]-yt)*(f[3]-yt));
            //Log.i("Debug", "distance: "+distance);
            if (distance<previous_distance){
                previous_distance = distance;
                result[0] = f[0];
                result[1] = f[1];
            }
        }
        return result;
    }
    @Override
    protected void onRestart() {
        super.onRestart();
        tabLayout.getTabAt(0).select();
    }

    public double[] rotateVector(double[] vector, double thetaDegrees) {
        double theta = Math.toRadians(thetaDegrees);  // Convert degrees to radians
        double cosTheta = Math.cos(theta);
        double sinTheta = Math.sin(theta);

        double xNew = cosTheta * vector[0] - sinTheta * vector[1];
        double yNew = sinTheta * vector[0] + cosTheta * vector[1];

        return new double[] {xNew, yNew};
    }

    public String pointDirection(double[] vector, double[] point) {
        double crossProduct = vector[0] * point[1] - vector[1] * point[0];

        if (crossProduct > 0) {
            Log.i("major", " left side");
            return "left";
        } else if (crossProduct < 0) {
            Log.i("major", " right side");
            return "right";
        } else {
            return "on the vector";
        }
    }

    public static int findMaxSpacing(List<double[]> points, int bufferSize) {
        double maxMeanDistance = 0;
        int maxIndex = 0;

        // Iterate over the points, keeping a sliding window of size 'bufferSize'
        for (int i = 0; i <= points.size() - bufferSize; i++) {
            double sumDistance = 0;

            // Calculate the sum of distances within the current window
            for (int j = i; j < i + bufferSize - 1; j++) {
                sumDistance += calculateDistance(points.get(j), points.get(j + 1));
            }

            // Calculate the mean distance in the window
            double meanDistance = sumDistance / (bufferSize - 1);

            // Update the max mean distance and index if needed
            if (meanDistance > maxMeanDistance) {
                maxMeanDistance = meanDistance;
                maxIndex = i;
            }
        }

        return maxIndex; // Return the index of the window with the maximum mean distance
    }
    private static double calculateDistance(double[] p1, double[] p2) {
        return Math.sqrt(Math.pow(p2[0] - p1[0], 2) + Math.pow(p2[1] - p1[1], 2));
    }

}