typedef struct
{
    double radius;  // Radius of ray in polar coordinates
    double angle;   // Angle of ray in polar coordinates
    double dRadius; // Rate of change of radius
} RayState;         // Rate of change of angle omitted - derived as `L/r^2`