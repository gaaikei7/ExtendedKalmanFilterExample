# Extended Kalman Filter Example
This is an simple implementation of the EKF over the bicycle model (aka Ackermann steering model).

Ackerman steering model is the next:

x/dt = v*cos(theta)
y/dt = v*sin(theta)
theta/dt = v/L*tan(phi)

where x is the position in the x-axis
y is the position in the y-axis
theta is the orientation of the vehicle
v is the speed of the vehicle
phi is the steering wheel angle
and L is the wheel base
