#include <math.h>
#include <pthread.h>
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

int n_thd;  // number of threads

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

void update_position(double* x,
                     double* y,
                     double* vx,
                     double* vy,
                     int n,
                     int thd,
                     pthread_mutex_t* mutexs) {
    // TODO: update position
    for (int i = thd * n; (i < (thd + 1) * n) && (i < n_body); ++i) {
        // printf("Thread %d Body %d: update_position\n", thd, i);
        pthread_mutex_lock(&mutexs[i]);
        if (x[i] < 0 || x[i] > bound_x) {
            vx[i] = -vx[i];
        }
        if (y[i] < 0 || y[i] > bound_y) {
            vy[i] = -vy[i];
        }
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
        pthread_mutex_unlock(&mutexs[i]);
    }
}

void update_velocity(double* m,
                     double* x,
                     double* y,
                     double* vx,
                     double* vy,
                     int n,
                     int thd,
                     pthread_mutex_t* mutexs) {
    // TODO: calculate force and acceleration, update velocity
    for (int i = thd * n; (i < (thd + 1) * n) && (i < n_body); ++i) {
        // printf("Thread %d Body %d: update_velocity\n", thd, i);
        double fx = 0, fy = 0;
        for (int j = 0; j < n_body; ++j) {
            if (i != j) {
                double ux = 0, uy = 0;
                pthread_mutex_lock(&mutexs[i]);
                double x0 = x[i], y0 = y[i];
                double vx0 = vx[i], vy0 = vy[i];
                pthread_mutex_unlock(&mutexs[i]);

                pthread_mutex_lock(&mutexs[j]);
                double x1 = x[j], y1 = y[j];
                double vx1 = vx[j], vy1 = vy[j];
                pthread_mutex_unlock(&mutexs[j]);

                double r2 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);

                if (r2 > radius2) {
                    ux = (x1 - x0) / sqrt(r2);
                    uy = (y1 - y0) / sqrt(r2);
                    double f = gravity_const * m[i] * m[j] / (r2 + err);
                    fx += f / m[i] * ux;
                    fy += f / m[i] * uy;

                } else if ((vx0 - vx1) * (x0 - x1) + (vy0 - vy1) * (y0 - y1) <
                           0) {
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
                }
                pthread_mutex_lock(&mutexs[i]);
                x[i] = x0, y[i] = y0;
                vx[i] = vx0, vy[i] = vy0;
                pthread_mutex_unlock(&mutexs[i]);

                pthread_mutex_lock(&mutexs[j]);
                x[j] = x1, y[j] = y1;
                vx[j] = vx1, vy[j] = vy1;
                pthread_mutex_unlock(&mutexs[j]);
            }
        }

        pthread_mutex_lock(&mutexs[i]);
        vx[i] += fx * dt;
        vy[i] += fy * dt;
        // printf("vx: %lf, vy: %lf\n", vx[i], vy[i]);
        pthread_mutex_unlock(&mutexs[i]);
    }
}

typedef struct {
    // TODO: specify your arguments for threads
    // int a;
    // int b;
    // TODO END
    int n = 0;
    double* m = nullptr;
    double* x = nullptr;
    double* y = nullptr;
    double* vx = nullptr;
    double* vy = nullptr;
    pthread_mutex_t* mutexs = nullptr;
    int thd = 0;
} Args;

void* worker(void* args) {
    // TODO: procedure in each threads

    Args* my_arg = (Args*)args;
    // int a = my_arg->a;
    // int b = my_arg->b;
    int thd = my_arg->thd;
    int n = my_arg->n;
    double* m = my_arg->m;
    double* x = my_arg->x;
    double* y = my_arg->y;
    double* vx = my_arg->vx;
    double* vy = my_arg->vy;
    pthread_mutex_t* mutexs = my_arg->mutexs;
    // printf("Thread %d: hello\n", thd);
    update_velocity(m, x, y, vx, vy, n, thd, mutexs);

    update_position(x, y, vx, vy, n, thd, mutexs);
    return NULL;
    // TODO END
}

Args* args;

void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];
    args = new Args[n_body];
    pthread_mutex_t* mutexs = new pthread_mutex_t[n_body];
    pthread_t* thds = new pthread_t[n_thd];
    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("sequential", n_body, bound_x, bound_y);
    for (int i = 0; i < n_body; i++) {
        pthread_mutex_init(&mutexs[i], NULL);
    }
    int block_size = (n_body - 1) / n_thd + 1;

    for (int i = 0; i < n_iteration; i++) {
        std::chrono::high_resolution_clock::time_point t1 =
            std::chrono::high_resolution_clock::now();
        // TODO: assign jobs

        for (int thd = 0; thd < n_thd; thd++) {
            args[thd].n = block_size;
            args[thd].thd = thd;
            args[thd].m = m;
            args[thd].x = x;
            args[thd].y = y;
            args[thd].vx = vx;
            args[thd].vy = vy;
            args[thd].mutexs = mutexs;
        }
        for (int thd = 0; thd < n_thd; thd++)
            pthread_create(&thds[thd], NULL, worker, &args[thd]);
        for (int thd = 0; thd < n_thd; thd++)
            pthread_join(thds[thd], NULL);
        // TODO End

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
    delete[] mutexs;
    delete[] args;
}

int main(int argc, char* argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);

#ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Pthread");
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(0, bound_x, 0, bound_y);
#endif
    master();

    return 0;
}
