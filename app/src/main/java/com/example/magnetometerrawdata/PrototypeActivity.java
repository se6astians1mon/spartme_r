package com.example.magnetometerrawdata;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Context;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.media.MediaPlayer;
import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.ImageView;
import android.widget.VideoView;

import com.google.android.material.tabs.TabLayout;
import com.opencsv.CSVReader;

import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;

public class PrototypeActivity extends AppCompatActivity {

    private TabLayout tabLayout;
    private Context context;
    private Button btn_calibrate;
    private ImageView img_map;
    private VideoView vid_map;
    private LinkedList<Float> samplesX;
    private LinkedList<Float> samplesY;
    private LinkedList<Float> samplesZ;
    private int bufferLength;
    float[] calibrationValues;
    float[] currPos;
    private SensorEventListener magnetSensorListener;
    private Sensor mag;
    private ArrayList<Float[]> mapping;
    private float offsetY;
    private float offsetX;
    private float scale;
    private boolean calibrated;


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.prototype);

        SensorManager sensorManager = (SensorManager) getSystemService(SENSOR_SERVICE);
        samplesX = new LinkedList<Float>();
        samplesY = new LinkedList<Float>();
        samplesZ = new LinkedList<Float>();
        bufferLength = 10;
        //values depending on smartphone screen and position of magnetometer
        scale = 140.0f;             //scale for step size movement
        offsetX = -1400;           //smartphone: left right
        offsetY = -2350;        //smartphone: down-up

        //will be set later by user
        calibrationValues = new float[]{0.0f, 0.0f, 0.0f};
        currPos = new float[]{0.0f,0.0f};
        //init queue
        for(int i = 0;i<bufferLength;i++){
            samplesX.add(0.0f);
            samplesY.add(0.0f);
            samplesZ.add(0.0f);
        }
        magnetSensorListener = new SensorEventListener() {
            @Override
            public void onSensorChanged(SensorEvent event) {
                if (event.sensor.getType() == Sensor.TYPE_MAGNETIC_FIELD) {
                    //calc running sum
                    samplesX.remove();
                    samplesY.remove();
                    samplesZ.remove();
                    samplesX.add(event.values[0]);
                    samplesY.add(event.values[1]);
                    samplesZ.add(event.values[2]);
                    //calculate position if it's calibrated :)
                    if(calibrated){
                        float [] pos = calcPos(mapping, avg(samplesX)-calibrationValues[0], avg(samplesY)-calibrationValues[1]);
                        if(pos[0]!=currPos[0] || pos[1]!=currPos[1]){
                            currPos[0]=pos[0];
                            currPos[1]=pos[1];
                            //update GUI
                            img_map.setY(currPos[1]*scale+offsetY);
                            img_map.setX(currPos[0]*scale+offsetX);
                            vid_map.setY(currPos[1]*scale+offsetY);
                            vid_map.setX(currPos[0]*scale+offsetX);
                        }   
                    }

                }
            }

            @Override
            public void onAccuracyChanged(Sensor sensor, int i) {

            }
        };
        mag = sensorManager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD_UNCALIBRATED);
        context = this;
        tabLayout = findViewById(R.id.tabLayout);
        tabLayout.getTabAt(1).select();
        tabLayout.addOnTabSelectedListener(new TabLayout.OnTabSelectedListener() {
            @Override
            public void onTabSelected(TabLayout.Tab tab) {
                Log.i("info", "chouette " + tab.getText());
                if(tab.getText().equals("development")){
                    Log.i("info", "chouette equals acc");
                    finish();
                }
            }

            @Override
            public void onTabUnselected(TabLayout.Tab tab) {

            }

            @Override
            public void onTabReselected(TabLayout.Tab tab) {

            }
        });

        img_map = findViewById(R.id.img_map);

        //video setup
        vid_map = findViewById(R.id.vid_map);
        vid_map.setScaleX(5.5f);
        vid_map.setScaleY(5.5f);
        String path = "android.resource://" + getPackageName() + "/" + R.raw.video;
        vid_map.setVideoPath(path);
        vid_map.setOnPreparedListener(new MediaPlayer.OnPreparedListener() {
            @Override
            public void onPrepared(MediaPlayer mp) {
                mp.setLooping(true);
            }
        });
        vid_map.start();

        btn_calibrate = findViewById(R.id.btn_calibrate);
        btn_calibrate.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                calibrationValues[0]= avg(samplesX);
                calibrationValues[1]= avg(samplesY);
                calibrationValues[2]= avg(samplesZ);
                btn_calibrate.setText("Calibrated");
                calibrated = true;
            }
        });
        //get mapping
        mapping = importMapping();
        //get it before we start measuring otherwise we use null object mapping.
        sensorManager.registerListener(magnetSensorListener, sensorManager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD), 1000)        ;
    }

    ArrayList<Float[]> importMapping(){
        //loads a mapping from real coordinates x [0] y [1] to magnetic field strength values x [2] y [3]
        ArrayList<Float[]> result = new ArrayList<>();
        try {
            Log.i("Info", "test");
            CSVReader reader = new CSVReader(new InputStreamReader(getResources().openRawResource(R.raw.magnetometer1cm20x33_minusearth_interpolated)));
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
                Log.i("info", "values: "+xValue + ", "+yValue+", "+xReal+", "+yReal);
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

    private float[] calcPos(ArrayList<Float[]> posList, float xt, float yt) {
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
}