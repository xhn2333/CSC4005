#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#ifdef GUI
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "./headers/logger.h"
#include "./headers/physics.h"

int n_body;
int n_iteration;

void generate_data(double* m,
                   double* x,
                   double* y,
                   double* vx,
                   double* vy,
                   int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}

void update_position(double* x, double* y, double* vx, double* vy, int n) {
    // TODO: update position
    for (int i = 0; i < n; ++i) {
        if (x[i] < 0 || x[i] > bound_x) {
            vx[i] = -vx[i];
        }
        if (y[i] < 0 || y[i] > bound_y) {
            vy[i] = -vy[i];
        }
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
    }
}

void update_velocity(double* m,
                     double* x,
                     double* y,
                     double* vx,
                     double* vy,
                     int n) {
    // TODO: calculate force and acceleration, update velocity

    for (int i = 0; i < n; ++i) {
        double fx = 0, fy = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double ux = 0, uy = 0;
                double r2 = (x[j] - x[i]) * (x[j] - x[i]) +
                            (y[j] - y[i]) * (y[j] - y[i]);
                if (r2 > radius2) {
                    ux = (x[j] - x[i]) / sqrt(r2);
                    uy = (y[j] - y[i]) / sqrt(r2);
                    double f = gravity_const * m[i] * m[j] / (r2 + err);
                    fx += f / m[i] * ux;
                    fy += f / m[i] * uy;

                } else if ((vx[i] - vx[j]) * (x[i] - x[j]) +
                               (vy[i] - vy[j]) * (y[i] - y[j]) <
                           0) {
                    double tmpvx = vx[i], tmpvy = vy[i];
                    vx[i] = ((m[i] - m[j]) * vx[i] + 2 * m[j] * vx[j]) /
                            (m[i] + m[j] + err);
                    vy[i] = ((m[i] - m[j]) * vy[i] + 2 * m[j] * vy[j]) /
                            (m[i] + m[j] + err);
                    vx[j] = (2 * m[i] * tmpvx + (m[j] - m[i]) * vx[j]) /
                            (m[i] + m[j] + err);
                    vy[j] = (2 * m[i] * tmpvy + (m[j] - m[i]) * vy[j]) /
                            (m[i] + m[j] + err);
                    // vx[i] = -vx[i];
                    // vy[i] = -vy[i];
                }
            }
        }
        vx[i] += fx * dt;
        vy[i] += fy * dt;
        // printf("vx: %lf, vy: %lf\n", vx[i], vy[i]);
    }
}

void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("sequential", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++) {
        std::chrono::high_resolution_clock::time_point t1 =
            std::chrono::high_resolution_clock::now();

        update_velocity(m, x, y, vx, vy, n_body);
        update_position(x, y, vx, vy, n_body);

        std::chrono::high_resolution_clock::time_point t2 =
            std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(x, y);

#ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++) {
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
#else

#endif
    }
    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
}

int main(int argc, char* argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

#ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
#endif
    master();

    printf("Student ID: 120090453\n");  // replace it with your student id
    printf("Name: Haonan XUE\n");       // replace it with your name
    printf("Assignment 3: N Body Simulation Sequential Implementation\n");

    return 0;
}
