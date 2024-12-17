# Magnetometer Prototype
A project aiming at identifying device position with centimeter precision. A magnet is used to create a magnetic field. The device's magnetometer is used to determine the device's position by magnetic fingerprinting.

## Code structure
All relevant code resides within the class "PrototypeActivity.java" in app/src/main/java/com/example/magnetometerrawdata. "MainActivity.java" allows for mapping a magnetic field from a custom magnet.

## Algorithm
The user initially places its device in parallel to the map and the magnet, but far away enough for the magnet not having any influence on the accelerometer. Clicking on "Calibrate" then uses the recorded
values to determine the orientation of the map with respect to the earth's magnetic field and subtract it whenever a measure is taken within the artificial magnetic field of the magnet.
The application then loads a mapping present in "app\src\main\res\raw\magnetometer1cm20x33_minusearth.csv". This file contains a mapping of the magnetic field of the magnet used. 
The magnetometer data is then smoothed by a moving average (line 75 and following). This value pair of magnetic field strength values is then looked up within the previously loaded mapping (line 182 and following). This is achieved by calculating the Euclidian distance for every point in the mapping to the value pair. The closest point is then used to move the map (minus the display, sensor and screen size offset, line 88 and 89)  

## Recording
The recording part of the application allows to record magnetometer values and save them with the corresponding position in a mysql table. An apache and mysql server are required. 
The mysql database named "situlearn" and a table "magnetometer" with the following columns (name type) must be created beforehand: 
(position_name 	text )  
(x 	            float )  
(y 	            float )  
(z 	            float )  
(realX 	        int(11))  
(realY 	        int(11))  
(mduration 	    int(11) )  
(minterval    	  int(11) )  
(time 	          bigint(20) 	)  

It is also recommended to print a numerated grid with the magnet on a fixed position. The device is then placed on this grid. Once the correct ip of the apache server serving db_access.php
(eg. 192.168.1.12/situlearn_db_access/db_access.php) is set alongside the number of values taken per seconds and the number of seconds, clicking "Record" records the values and sends them to
the database. Once the mapping is complete, one can export it with a GROUP BY sql query in order to obtain unique magnetic field values for each position. A last measure has to be taken where
the device is placed in parallel outside the magnet's field. This value has then to be subtracted from the data obtained. The result is the field of the magnet without the influence of the
earth's magnetic field.


## Limitations
The use is limited to the reach of the magnet. Using a strong magnet can lead to magnetisation of parts of the device itself, temporarily falsifying all measurements. The precision decreases over 
distance: While the gap between magnetic field values close to the magnet is large, the farer away from the magnet, the closer the values are next to each other, reducing precision. The loss of 
precision decreases cubically with distance.

## Mapping
More precise mappings could potentially obtained from simulations, using software such as FEMM. The magnet which has been used in this experiment can be found in "scripts_and_models".
To export the data, modify and execute exportFemm.lua from within FEMM: x1,x2,y1,y2 are two points of a rectangle defining the area to be exported. dx and dy define the distance between
any two points (eg. 0,1 equals 1mm precision) and line 19 defines the output file.
