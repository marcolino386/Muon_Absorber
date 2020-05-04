set terminal png size 500,500
set parametric
set trange[0:2*pi]
set size square
set output 'hits_detector.png'

z_0 = 0.9
carbon_pDz = 1.125
concrete_pDz = 0.76
pos_after_detec = 3
initial_angle = 2
final_angle = 10

initial_radius = (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec)*tan(initial_angle*pi/180.00)
final_radius = (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec)*tan(final_angle*pi/180.00)

plot initial_radius*cos(t), initial_radius*sin(t) title " ", final_radius*cos(t) , final_radius*sin(t) title " ", "position_mu_p.csv" u 1:2 title "mu+", "position_mu_m.csv" u 1:2 title "mu-"

