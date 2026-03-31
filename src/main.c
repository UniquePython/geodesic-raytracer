#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#include <raylib.h>

#define WIDTH 900  // Window width
#define HEIGHT 900 // Window height

#include <pthread.h>

#define NUM_THREADS 20 // My machine has 20 cores

typedef struct
{
    int startRow;
    int endRow;
    Color *pixelBuffer;
    Image *starfield;
} ThreadData;

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
    ESCAPED,  // Ray escaped black hole
    HIT_DISK, // Hit accretion disk
} Outcome;

typedef struct trace_result
{
    Outcome outcome;
    RayState finalState;
} TraceResult;

#define RS 1.0     // Schwarzschild radius, normalized to 1
#define R_CAM 30.0 // Camera distance from black hole
#define R_MAX 50.0 // Stands for infinity

#define FOV (PI / 3.0)       // 60 deg
#define CAM_INCLINATION 30.0 // 30 deg

#define DLAMBDA 0.1
#define MAX_STEPS 10000

RayState derivatives(RayState currState, double angularMomentum, double carterConst)
{
    RayState derivative;

    derivative.radius = currState.dRadius;

    double radiusSquared = currState.radius * currState.radius;
    double angularMomentumSquared = angularMomentum * angularMomentum;

    derivative.angle = angularMomentum / radiusSquared;

    double effectiveLSq = angularMomentumSquared + carterConst;
    double outwardCentrifugalPush = effectiveLSq / (radiusSquared * currState.radius);
    double inwardPull = 3 * effectiveLSq * RS / (2 * radiusSquared * radiusSquared);

    derivative.dRadius = outwardCentrifugalPush - inwardPull;

    derivative.theta = currState.dTheta;

    double dragFromRadialMotion = -(2.0 / currState.radius) * currState.dRadius * currState.dTheta;

    double sinTheta = sin(currState.theta);
    double cosTheta = cos(currState.theta);
    double couplingBwAzimuthalAndPolarMotion = sinTheta * cosTheta * (angularMomentum / radiusSquared) * (angularMomentum / radiusSquared);

    derivative.dTheta = dragFromRadialMotion + couplingBwAzimuthalAndPolarMotion;

    return derivative;
}

RayState rk4_step(RayState currState, double angularMomentum, double carterConst, double dLambda)
{
    RayState k1 = derivatives(currState, angularMomentum, carterConst);

    // --- K2 ------------>

    RayState k2Input = {
        .radius = currState.radius + (k1.radius * dLambda / 2),
        .dRadius = currState.dRadius + (k1.dRadius * dLambda / 2),
        .angle = currState.angle + (k1.angle * dLambda / 2),
        .theta = currState.theta + (k1.theta * dLambda / 2),
        .dTheta = currState.dTheta + (k1.dTheta * dLambda / 2),
    };

    RayState k2 = derivatives(k2Input, angularMomentum, carterConst);

    // --- K3 ------------>

    RayState k3Input = {
        .radius = currState.radius + (k2.radius * dLambda / 2),
        .dRadius = currState.dRadius + (k2.dRadius * dLambda / 2),
        .angle = currState.angle + (k2.angle * dLambda / 2),
        .theta = currState.theta + (k2.theta * dLambda / 2),
        .dTheta = currState.dTheta + (k2.dTheta * dLambda / 2),
    };

    RayState k3 = derivatives(k3Input, angularMomentum, carterConst);

    // --- K4 ------------>

    RayState k4Input = {
        .radius = currState.radius + (k3.radius * dLambda),
        .dRadius = currState.dRadius + (k3.dRadius * dLambda),
        .angle = currState.angle + (k3.angle * dLambda),
        .theta = currState.theta + (k3.theta * dLambda),
        .dTheta = currState.dTheta + (k3.dTheta * dLambda),
    };

    RayState k4 = derivatives(k4Input, angularMomentum, carterConst);

    // --- Final Blend ------------>

    RayState result = {
        .radius = currState.radius + (dLambda / 6) * (k1.radius + 2 * k2.radius + 2 * k3.radius + k4.radius),
        .angle = currState.angle + (dLambda / 6) * (k1.angle + 2 * k2.angle + 2 * k3.angle + k4.angle),
        .dRadius = currState.dRadius + (dLambda / 6) * (k1.dRadius + 2 * k2.dRadius + 2 * k3.dRadius + k4.dRadius),
        .theta = currState.theta + (dLambda / 6) * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta),
        .dTheta = currState.dTheta + (dLambda / 6) * (k1.dTheta + 2 * k2.dTheta + 2 * k3.dTheta + k4.dTheta),
    };

    return result;
}

TraceResult trace_ray(RayState initial, double angularMomentum, double carterConst, double dLambda, int maxSteps)
{
    RayState state = initial;
    RayState prevState = initial;

    for (int stepsTaken = 0; stepsTaken < maxSteps; stepsTaken++)
    {
        prevState = state;

        double distToHorizon = state.radius - RS;
        double distToPhotonSphere = fabs(state.radius - 1.5 * RS);
        double minDist = fmin(distToHorizon, distToPhotonSphere);

        double adaptiveLambda = dLambda * fmin(1.0, minDist / (3.0 * RS));
        adaptiveLambda = fmax(0.05, adaptiveLambda);

        state = rk4_step(state, angularMomentum, carterConst, adaptiveLambda);

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

        bool crossedEquatorialPlane = (prevState.theta - PI / 2) * (state.theta - PI / 2) < 0;
        bool inDiskRegion = state.radius >= 3.0 * RS && state.radius <= 15.0 * RS;

        if (crossedEquatorialPlane && inDiskRegion)
        {
            double t = (PI / 2.0 - prevState.theta) / (state.theta - prevState.theta);

            RayState interpolated = {
                .radius = prevState.radius + t * (state.radius - prevState.radius),
                .angle = prevState.angle + t * (state.angle - prevState.angle),
                .theta = PI / 2.0,
                .dRadius = prevState.dRadius + t * (state.dRadius - prevState.dRadius),
                .dTheta = prevState.dTheta + t * (state.dTheta - prevState.dTheta),
            };

            TraceResult result = {
                .outcome = HIT_DISK,
                .finalState = interpolated,
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

void *render_rows(void *arg)
{
    ThreadData *data = (ThreadData *)arg;

    for (int y = data->startRow; y < data->endRow; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            double screenX = (x - WIDTH / 2.0) / WIDTH;
            double screenY = (y - HEIGHT / 2.0) / HEIGHT;

            double ax = screenX * FOV;
            double ay = screenY * FOV;
            double angularMomentum = R_CAM * sin(ax);

            double cameraInclination = CAM_INCLINATION * (PI / 180.0);

            RayState initial = {
                .radius = R_CAM,
                .angle = 0.0,
                .theta = (PI / 2.0) - cameraInclination,
                .dRadius = -1.0,
                .dTheta = sin(ay) / R_CAM,
            };

            double sinTheta = sin(initial.theta);
            double cosTheta = cos(initial.theta);

            double carterConst;
            if (fabs(sinTheta) < 1e-6)
                carterConst = R_CAM * R_CAM * R_CAM * R_CAM * (initial.dTheta * initial.dTheta);
            else
                carterConst = R_CAM * R_CAM * R_CAM * R_CAM * (initial.dTheta * initial.dTheta) + (cosTheta * cosTheta) / (sinTheta * sinTheta) * angularMomentum * angularMomentum;

            TraceResult result = trace_ray(initial, angularMomentum, carterConst, DLAMBDA, MAX_STEPS);

            Color color;
            if (result.outcome == CAPTURED)
            {
                color = BLACK;
            }
            else if (result.outcome == HIT_DISK)
            {
                double r = result.finalState.radius;
                double innerEdge = 3.0 * RS;
                double outerEdge = 15.0 * RS;

                double t = (r - innerEdge) / (outerEdge - innerEdge);
                t = fmax(0.0, fmin(1.0, t));

                Color diskColor;
                if (t < 0.25)
                {
                    double s = t / 0.25;
                    diskColor = (Color){255, (unsigned char)(255 - s * 55), (unsigned char)(255 - s * 255), 255};
                }
                else if (t < 0.6)
                {
                    double s = (t - 0.25) / 0.35;
                    diskColor = (Color){255, (unsigned char)(200 - s * 130), 0, 255};
                }
                else
                {
                    double s = (t - 0.6) / 0.4;
                    diskColor = (Color){(unsigned char)(255 - s * 155), (unsigned char)(70 - s * 70), 0, 255};
                }

                double brightness = 1.0 / sqrt(r);
                brightness = fmax(0.1, fmin(1.0, brightness));
                diskColor.r = (unsigned char)(diskColor.r * brightness);
                diskColor.g = (unsigned char)(diskColor.g * brightness);
                diskColor.b = (unsigned char)(diskColor.b * brightness);

                color = diskColor;
            }
            else
            {
                double u = fmod(result.finalState.angle, 2 * PI) / (2 * PI);
                if (u < 0)
                    u += 1.0;

                double v = fmax(0.0, fmin(1.0, result.finalState.theta / PI));

                int texX = (int)(u * data->starfield->width);
                int texY = (int)(v * data->starfield->height);

                texX = (int)fmax(0, fmin(texX, data->starfield->width - 1));
                texY = (int)fmax(0, fmin(texY, data->starfield->height - 1));

                color = GetImageColor(*data->starfield, texX, texY);
            }

            data->pixelBuffer[(HEIGHT - 1 - y) * WIDTH + x] = color;
        }
    }

    return NULL;
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

    static Color pixelBuffer[WIDTH * HEIGHT];

    bool rendered = false;

    while (!WindowShouldClose())
    {
        if (!rendered)
        {
            pthread_t threads[NUM_THREADS];
            ThreadData threadData[NUM_THREADS];

            int rowsPerThread = HEIGHT / NUM_THREADS;

            for (int i = 0; i < NUM_THREADS; i++)
            {
                threadData[i].startRow = i * rowsPerThread;
                threadData[i].endRow = (i == NUM_THREADS - 1) ? HEIGHT : (i + 1) * rowsPerThread;
                threadData[i].starfield = &starfield;
                threadData[i].pixelBuffer = pixelBuffer;
                pthread_create(&threads[i], NULL, render_rows, &threadData[i]);
            }

            for (int i = 0; i < NUM_THREADS; i++)
                pthread_join(threads[i], NULL);

            UpdateTexture(target.texture, pixelBuffer);

            rendered = true;
        }

        BeginDrawing();
        DrawTextureRec(target.texture, (Rectangle){0, 0, WIDTH, -HEIGHT}, (Vector2){0, 0}, WHITE);
        EndDrawing();
    }

    CloseWindow();
    UnloadImage(starfield);

    return 0;
}