set cp=.;jfreechart-1.0.13.jar;itext-1.4.3.jar;jcommon-1.0.16.jar
set infiles=data/Embryo.txt data/L1.txt data/L2.txt data/L3.txt data/L4.txt data/Young_Adult_Male.txt data/Young_Adult.txt data/Day0.txt data/Day5.txt data/Day8.txt data/Day12.txt

java -Xmx256m -classpath %cp% SaturationPlotsCreator output_all_abs.pdf output_all_abs.txt 0        0   %infiles%
java -Xmx256m -classpath %cp% SaturationPlotsCreator output_all_fra.pdf output_all_fra.txt 61569834 0   %infiles%
java -Xmx256m -classpath %cp% SaturationPlotsCreator output_sam_abs.pdf output_sam_abs.txt 0        100 %infiles%
java -Xmx256m -classpath %cp% SaturationPlotsCreator output_sam_fra.pdf output_sam_fra.txt 61569834 100 %infiles%
