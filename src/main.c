typedef struct
{
    double radius;  // Radius of ray in polar coordinates
    double angle;   // Angle of ray in polar coordinates
    double dRadius; // Rate of change of radius
} RayState;         // Rate of change of angle omitted - derived as `L/r^2`

#define RS 1.0     // Schwarzschild radius, normalized to 1
#define R_CAM 20.0 // camera distance from black hole