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

RayState rk4_step(RayState currState, double angularMomentum, double dlambda)
{
    RayState k1 = derivatives(currState, angularMomentum);

    // --- K2 ------------>

    RayState k2Input = {
        .radius = currState.radius + (k1.radius * dlambda / 2),
        .dRadius = currState.dRadius + (k1.dRadius * dlambda / 2),
        .angle = currState.angle + (k1.angle * dlambda / 2),
    };

    RayState k2 = derivatives(k2Input, angularMomentum);

    // --- K3 ------------>

    RayState k3Input = {
        .radius = currState.radius + (k2.radius * dlambda / 2),
        .dRadius = currState.dRadius + (k2.dRadius * dlambda / 2),
        .angle = currState.angle + (k2.angle * dlambda / 2),
    };

    RayState k3 = derivatives(k3Input, angularMomentum);

    // --- K4 ------------>

    RayState k4Input = {
        .radius = currState.radius + (k3.radius * dlambda),
        .dRadius = currState.dRadius + (k3.dRadius * dlambda),
        .angle = currState.angle + (k3.angle * dlambda),
    };

    RayState k4 = derivatives(k4Input, angularMomentum);

    // --- Final Blend ------------>

    RayState result = {
        .radius = currState.radius + (dlambda / 6) * (k1.radius + 2 * k2.radius + 2 * k3.radius + k4.radius),
        .angle = currState.angle + (dlambda / 6) * (k1.angle + 2 * k2.angle + 2 * k3.angle + k4.angle),
        .dRadius = currState.dRadius + (dlambda / 6) * (k1.dRadius + 2 * k2.dRadius + 2 * k3.dRadius + k4.dRadius),
    };

    return result;
}