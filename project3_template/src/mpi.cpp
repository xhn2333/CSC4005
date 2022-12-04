#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <new>
#include <utility>
#include "./headers/logger.h"
#include "./headers/physics.h"

#ifdef GUI
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

int n_body;
int n_iteration;

int my_rank;
int world_size;

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
    for (int i = my_rank * n; (i < (my_rank + 1) * n) && (i < n_body); ++i) {
        // printf("Thread %d Body %d: update_position\n", thd, i);

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
    for (int i = my_rank * n; (i < (my_rank + 1) * n) && (i < n_body); ++i) {
        // printf("Thread %d Body %d: update_velocity\n", thd, i);
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
                    */
                    vx[i] = -vx[i];
                    vy[i] = -vy[i];
                }

                x[i] = x0, y[i] = y0;
                vx[i] = vx0, vy[i] = vy0;

                x[j] = x1, y[j] = y1;
                vx[j] = vx1, vy[j] = vy1;
            }
        }

        vx[i] += fx * dt;
        vy[i] += fy * dt;
        // printf("vx: %lf, vy: %lf\n", vx[i], vy[i]);
    }
}

int main(int argc, char* argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int block_size = (n_body - 1) / world_size + 1;

    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];
    int sendcount[n_body];
    int displs[n_body];

    if (my_rank == 0) {
#ifdef GUI
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
        glutInitWindowSize(500, 500);
        glutInitWindowPosition(0, 0);
        glutCreateWindow("N Body Simulation MPI Implementation");
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glMatrixMode(GL_PROJECTION);
        gluOrtho2D(0, bound_x, 0, bound_y);
#endif
    }
    generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);
    Logger* l = nullptr;
    if (my_rank == 0)
        l = new Logger("sequential", n_body, bound_x, bound_y);

    for (int i = 0; i < world_size; ++i) {
        displs[i] = i * block_size;
        sendcount[i] =
            displs[i] + block_size > n_body ? n_body - displs[i] : block_size;
    }

    for (int i = 0; i < n_iteration; i++) {
        std::chrono::high_resolution_clock::time_point t1, t2;
        if (my_rank == 0) {
            t1 = std::chrono::high_resolution_clock::now();
        }
        double tmpm[n_body];
        double tmpx[n_body];
        double tmpy[n_body];
        double tmpvx[n_body];
        double tmpvy[n_body];
        MPI_Allgatherv(total_m + displs[my_rank], sendcount[my_rank],
                       MPI_DOUBLE, tmpm, sendcount, displs, MPI_DOUBLE,
                       MPI_COMM_WORLD);
        MPI_Allgatherv(total_x + displs[my_rank], sendcount[my_rank],
                       MPI_DOUBLE, tmpx, sendcount, displs, MPI_DOUBLE,
                       MPI_COMM_WORLD);
        MPI_Allgatherv(total_y + displs[my_rank], sendcount[my_rank],
                       MPI_DOUBLE, tmpy, sendcount, displs, MPI_DOUBLE,
                       MPI_COMM_WORLD);
        MPI_Allgatherv(total_vx + displs[my_rank], sendcount[my_rank],
                       MPI_DOUBLE, tmpvx, sendcount, displs, MPI_DOUBLE,
                       MPI_COMM_WORLD);
        MPI_Allgatherv(total_vy + displs[my_rank], sendcount[my_rank],
                       MPI_DOUBLE, tmpvy, sendcount, displs, MPI_DOUBLE,
                       MPI_COMM_WORLD);

        for (int i = 0; i < n_body; ++i) {
            total_m[i] = tmpm[i];
            total_x[i] = tmpx[i];
            total_y[i] = tmpy[i];
            total_vx[i] = tmpvx[i];
            total_vy[i] = tmpvy[i];
        }

        // TODO: MPI routine
        update_velocity(total_m, total_x, total_y, total_vx, total_vy,
                        block_size);
        update_position(total_x, total_y, total_vx, total_vy, block_size);

        // TODO End
        if (my_rank == 0) {
            t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = t2 - t1;
            printf("Iteration %d, elapsed time: %.3f\n", i, time_span);
            (*l).save_frame(total_x, total_y);

#ifdef GUI
            glClear(GL_COLOR_BUFFER_BIT);
            glColor3f(1.0f, 0.0f, 0.0f);
            glPointSize(2.0f);
            glBegin(GL_POINTS);
            double xi;
            double yi;
            for (int i = 0; i < n_body; i++) {
                xi = total_x[i];
                yi = total_y[i];
                glVertex2f(xi, yi);
            }
            glEnd();
            glFlush();
            glutSwapBuffers();
#else

#endif
        }
    }

    delete[] total_m;
    delete[] total_x;
    delete[] total_y;
    delete[] total_vx;
    delete[] total_vy;
    if (my_rank == 0) {
        printf("Student ID: 120090453\n");  // replace it with your student id
        printf("Name: Haonan XUE\n");       // replace it with your name
        printf("Assignment 2: N Body Simulation MPI Implementation\n");
    }

    MPI_Finalize();

    return 0;
}
