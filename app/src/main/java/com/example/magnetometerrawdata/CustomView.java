package com.example.magnetometerrawdata;

import android.content.Context;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.util.AttributeSet;
import android.util.Log;
import android.view.View;

import java.util.Queue;

public class CustomView extends View {
    public Queue<float[]> bufferX;
    public int bufferLength;
    public double[] vector = {0.0, 0.0};
    public double[] center = {0.0, 0.0};
    float zoomfactor = 10;
    float length;
    float data[] = {0.0f, 0.0f, 0.0f};
    public float[] pointSpace = {0.0f, 0.0f};

    public CustomView(Context context, AttributeSet attr) {
        super(context, attr);
    }

    @Override
    public void onDraw(Canvas canvas) {
        zoomfactor = 300.0f/length;
        Log.i("zoom", zoomfactor+"");
        canvas.drawColor(Color.rgb(32, 32, 32));

        Paint redPaint = new Paint();
        redPaint.setColor(Color.rgb(255, 0, 0));

        Paint greenPaint = new Paint();
        greenPaint.setColor(Color.rgb(0, 255, 0));
        greenPaint.setStrokeWidth(4.0f);

        Paint whitePaint = new Paint();
        whitePaint.setColor(Color.rgb(255, 255, 255));

        Paint yellowPaint = new Paint();
        yellowPaint.setColor(Color.rgb(255, 255, 0));

        Paint pinkPaint = new Paint();
        pinkPaint.setColor(Color.rgb(255, 0, 255));

        Paint bluePaint = new Paint();
        bluePaint.setColor(Color.rgb(51, 153, 255));
        bluePaint.setStrokeWidth(4.0f);


        //draw cross hair
        canvas.drawLine(canvas.getWidth() / 2 + 0.0f, 0.0f, canvas.getWidth() / 2 + 0.0f, canvas.getHeight() + 0.0f, redPaint);
        canvas.drawLine(0.0f, canvas.getHeight() / 2 + 0.0f, canvas.getWidth() + 0.0f, canvas.getHeight() / 2 + 0.0f, redPaint);

        //draw vector lines
        //x over z
        //canvas.drawLine(canvas.getWidth() / 2, canvas.getHeight() / 2, canvas.getWidth() / 2 + data[0]*zoomfactor, canvas.getHeight() / 2 + data[2]*zoomfactor, pinkPaint);
        //x over y
        //draw fast side
        canvas.drawCircle(canvas.getWidth() / 2 - pointSpace[0]*zoomfactor,canvas.getHeight() / 2 + pointSpace[1]*zoomfactor, 15, redPaint);
        //draw center of eclipse
        float centerX = canvas.getWidth() / 2 - (float)center[0]*zoomfactor;
        float centerY = canvas.getHeight() / 2 + (float)center[1]*zoomfactor;
        canvas.drawCircle(centerX, centerY, 15, whitePaint);

        canvas.drawLine(centerX, centerY, centerX - ((float)vector[0])*400, centerY + ((float)vector[1])*400, yellowPaint);
        canvas.drawLine(canvas.getWidth() / 2, canvas.getHeight() / 2, canvas.getWidth() / 2 - data[0]*zoomfactor, canvas.getHeight() / 2 + data[1]*zoomfactor, greenPaint);
        // draw history
        for(float[] number : bufferX)
            canvas.drawCircle(canvas.getWidth() / 2 - number[0]*zoomfactor,canvas.getHeight() / 2 + number[1]*zoomfactor, 5, pinkPaint);
    }
}
