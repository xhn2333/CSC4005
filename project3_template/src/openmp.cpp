#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#ifdef GUI
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "./headers/logger.h"
#include "./headers/physics.h"

int n_body;
int n_iteration;

int n_omp_threads;

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

void update_position(double* x, double* y, double* vx, double* vy, int i) {
    // TODO: update position
    if (i < n_body) {
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
                     int i) {
    // TODO: calculate force and acceleration, update velocity
    if (i < n_body) {
        double fx = 0, fy = 0;
        for (int j = 0; j < n_body; ++j) {
            if (i != j) {
                double ux = 0, uy = 0;

                double x0 = x[i], y0 = y[i];
                double vx0 = vx[i], vy0 = vy[i];

                double x1 = x[j], y1 = y[j];
                double vx1 = vx[j], vy1 = vy[j];

                double r2 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);

                if (r2 > radius2) {
                    ux = (x1 - x0) / sqrt(r2);
                    uy = (y1 - y0) / sqrt(r2);
                    double f = gravity_const * m[i] * m[j] / (r2 + err);
                    fx += f / m[i] * ux;
                    fy += f / m[i] * uy;

                } else if ((vx0 - vx1) * (x0 - x1) + (vy0 - vy1) * (y0 - y1) <
                           0) {
                    /*
                    double tmpvx = vx0, tmpvy = vy0;
                    vx0 = ((m[i] - m[j]) * vx0 + 2 * m[j] * vx1) /
                          (m[i] + m[j] + err);
                    vy0 = ((m[i] - m[j]) * vy0 + 2 * m[j] * vy1) /
                          (m[i] + m[j] + err);
                    vx1 = (2 * m[i] * tmpvx + (m[j] - m[i]) * vx1) /
                          (m[i] + m[j] + err);
                    vy1 = (2 * m[i] * tmpvy + (m[j] - m[i]) * vy1) /
                          (m[i] + m[j] + err);
                    // vx[i] = -vx[i];
                    // vy[i] = -vy[i];
                    */
                    x[i] = -x0, y[i] = -y0;
                    vx[i] = -vx0, vy[i] = -vy0;
                }

                /*
                x[j] = x1, y[j] = y1;
                vx[j] = vx1, vy[j] = vy1;*/
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

        // TODO: choose better threads configuration
        omp_set_num_threads(n_omp_threads);
#pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_velocity(m, x, y, vx, vy, i);
        }

        omp_set_num_threads(n_omp_threads);
#pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_position(x, y, vx, vy, i);
        }

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
    n_omp_threads = atoi(argv[3]);

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

    printf("Student ID: 119010001\n");  // replace it with your student id
    printf("Name: Your Name\n");        // replace it with your name
    printf("Assignment 2: N Body Simulation OpenMP Implementation\n");

    return 0;
}
