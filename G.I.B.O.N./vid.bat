cd img
 ffmpeg -f image2 -framerate 90  -start_number 1 -i gibbon_%04d.png -s 1834x1051  -b:v 10000k ../res_1.avi  
cd ..
