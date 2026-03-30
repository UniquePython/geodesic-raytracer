#include <stdbool.h>
#include <stddef.h>

typedef struct
{
    double radius;  // Radius of ray in polar coordinates
    double angle;   // Angle of ray in polar coordinates
    double dRadius; // Rate of change of radius
} RayState;         // Rate of change of angle omitted - derived as `L/r^2` where L is angular momentum

typedef enum
{
    CAPTURED, // Ray fell inside black hole
    ESCAPED   // Ray escaped black hole
} Outcome;

#define RS 1.0     // Schwarzschild radius, normalized to 1
#define R_CAM 20.0 // Camera distance from black hole
#define R_MAX 50.0 // Stands for infinity

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

RayState rk4_step(RayState currState, double angularMomentum, double dLambda)
{
    RayState k1 = derivatives(currState, angularMomentum);

    // --- K2 ------------>

    RayState k2Input = {
        .radius = currState.radius + (k1.radius * dLambda / 2),
        .dRadius = currState.dRadius + (k1.dRadius * dLambda / 2),
        .angle = currState.angle + (k1.angle * dLambda / 2),
    };

    RayState k2 = derivatives(k2Input, angularMomentum);

    // --- K3 ------------>

    RayState k3Input = {
        .radius = currState.radius + (k2.radius * dLambda / 2),
        .dRadius = currState.dRadius + (k2.dRadius * dLambda / 2),
        .angle = currState.angle + (k2.angle * dLambda / 2),
    };

    RayState k3 = derivatives(k3Input, angularMomentum);

    // --- K4 ------------>

    RayState k4Input = {
        .radius = currState.radius + (k3.radius * dLambda),
        .dRadius = currState.dRadius + (k3.dRadius * dLambda),
        .angle = currState.angle + (k3.angle * dLambda),
    };

    RayState k4 = derivatives(k4Input, angularMomentum);

    // --- Final Blend ------------>

    RayState result = {
        .radius = currState.radius + (dLambda / 6) * (k1.radius + 2 * k2.radius + 2 * k3.radius + k4.radius),
        .angle = currState.angle + (dLambda / 6) * (k1.angle + 2 * k2.angle + 2 * k3.angle + k4.angle),
        .dRadius = currState.dRadius + (dLambda / 6) * (k1.dRadius + 2 * k2.dRadius + 2 * k3.dRadius + k4.dRadius),
    };

    return result;
}

Outcome trace_ray(RayState initial, double angularMomentum, double dLambda, int maxSteps)
{
    RayState state = initial;

    for (int stepsTaken = 0; stepsTaken < maxSteps; stepsTaken++)
    {
        state = rk4_step(state, angularMomentum, dLambda);

        bool crossedEventHorizon = state.radius <= RS;
        bool escapedToInfinity = state.radius > R_MAX;

        if (crossedEventHorizon)
            return CAPTURED;
        if (escapedToInfinity)
            return ESCAPED;
    }

    return ESCAPED;
}