set cbrange [0:18]
set palette rgbformulae 33,13,10
plot "vector_field.dat" using 1:2:3:4:5 with vectors linecolor palette z
