typedef struct
{
    double radius;  // Radius of ray in polar coordinates
    double angle;   // Angle of ray in polar coordinates
    double dRadius; // Rate of change of radius
} RayState;         // Rate of change of angle omitted - derived as `L/r^2` where L is angular momentum

#define RS 1.0     // Schwarzschild radius, normalized to 1
#define R_CAM 20.0 // camera distance from black hole

RayState derivatives(RayState currState, double angularMomentum)
{
    RayState derivative;

    derivative.radius = currState.dRadius;

    double radiusSquared = currState.radius * currState.radius;
    double angularMomentumSquared = angularMomentum * angularMomentum;

    derivative.angle = angularMomentum / radiusSquared;

    double outwardCentrifugalPush = angularMomentumSquared / (radiusSquared * currState.radius);
    double inwardPull = 3 * angularMomentumSquared * RS / (2 * radiusSquared * radiusSquared);

    derivative.dRadius = outwardCentrifugalPush - inwardPull;

    return derivative;
}