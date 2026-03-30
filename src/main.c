#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#include <raylib.h>

#define WIDTH 900  // Window width
#define HEIGHT 900 // Window height

typedef struct ray_state
{
    double radius;  // Radius of ray in polar coordinates
    double angle;   // Angle of ray in polar coordinates
    double theta;   // polar angle from north pole. PI/2 = equatorial plane
    double dRadius; // Rate of change of radius
    double dTheta;  // rate of change of theta
} RayState;         // Rate of change of angle omitted - derived as `L/r^2` where L is angular momentum

typedef enum outcome
{
    CAPTURED, // Ray fell inside black hole
    ESCAPED   // Ray escaped black hole
} Outcome;

typedef struct trace_result
{
    Outcome outcome;
    RayState finalState;
} TraceResult;

#define RS 1.0     // Schwarzschild radius, normalized to 1
#define R_CAM 20.0 // Camera distance from black hole
#define R_MAX 50.0 // Stands for infinity

#define FOV (PI / 3.0)       // 60 deg
#define CAM_INCLINATION 20.0 // 20 deg

#define DLAMBDA 0.1
#define MAX_STEPS 10000

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

    derivative.theta = currState.dTheta;

    double dragFromRadialMotion = -(2.0 / currState.radius) * currState.dRadius * currState.dTheta;

    double sinTheta = sin(currState.theta);
    double cosTheta = cos(currState.theta);
    double couplingBwAzimuthalAndPolarMotion = sinTheta * cosTheta * (angularMomentum / radiusSquared) * (angularMomentum / radiusSquared);

    derivative.dTheta = dragFromRadialMotion + couplingBwAzimuthalAndPolarMotion;

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

TraceResult trace_ray(RayState initial, double angularMomentum, double dLambda, int maxSteps)
{
    RayState state = initial;

    for (int stepsTaken = 0; stepsTaken < maxSteps; stepsTaken++)
    {
        state = rk4_step(state, angularMomentum, dLambda);

        bool crossedEventHorizon = state.radius <= RS;
        bool escapedToInfinity = state.radius > R_MAX;

        if (crossedEventHorizon)
        {
            TraceResult result = {
                .outcome = CAPTURED,
                .finalState = state,
            };

            return result;
        }
        if (escapedToInfinity)
        {
            TraceResult result = {
                .outcome = ESCAPED,
                .finalState = state,
            };

            return result;
        }
    }

    TraceResult result = {
        .outcome = ESCAPED,
        .finalState = state,
    };

    return result;
}

int main(void)
{
    InitWindow(WIDTH, HEIGHT, "Geodesic Ray Tracing in Curved Spacetime");
    SetTargetFPS(GetMonitorRefreshRate(GetCurrentMonitor()));

    RenderTexture2D target = LoadRenderTexture(WIDTH, HEIGHT);

    Image starfield = LoadImage("data/starmap_4k.jpg");
    if (starfield.data == NULL)
    {
        TraceLog(LOG_ERROR, "Failed to load starfield image");
        return 1;
    }

    int currentRow = 0;

    while (!WindowShouldClose())
    {
        if (currentRow < HEIGHT)
        {
            BeginTextureMode(target);
            for (int x = 0; x < WIDTH; x++)
            {
                double screenX = (x - WIDTH / 2.0) / WIDTH;
                double screenY = (currentRow - HEIGHT / 2.0) / HEIGHT;
                double alpha = sqrt(screenX * screenX + screenY * screenY) * FOV;
                double angularMomentum = R_CAM * sin(alpha);

                double cameraInclination = CAM_INCLINATION * (PI / 180.0);

                RayState initial = {
                    .radius = R_CAM,
                    .angle = 0.0,
                    .theta = (PI / 2.0) - cameraInclination + (screenY * FOV),
                    .dRadius = -1.0,
                    .dTheta = 0.0,
                };

                TraceResult result = trace_ray(initial, angularMomentum, DLAMBDA, MAX_STEPS);

                Color color;
                if (result.outcome == CAPTURED)
                    color = BLACK;
                else
                {
                    double u = fmod(result.finalState.angle, 2 * PI) / (2 * PI);
                    if (u < 0)
                        u += 1.0;
                    double v = currentRow / (double)HEIGHT;

                    int texX = (int)(u * starfield.width);
                    int texY = (int)(v * starfield.height);
                    color = GetImageColor(starfield, texX, texY);
                }

                DrawPixel(x, currentRow, color);
            }
            EndTextureMode();
            currentRow++;
        }

        BeginDrawing();
        DrawTextureRec(target.texture, (Rectangle){0, 0, WIDTH, -HEIGHT}, (Vector2){0, 0}, WHITE);
        EndDrawing();
    }

    CloseWindow();
    UnloadImage(starfield);

    return 0;
}