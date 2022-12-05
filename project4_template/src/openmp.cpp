#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <omp.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"

int size;
int n_omp_threads;
int n_iter;

void initialize(float* data, int i) {
    // intialize the temperature distribution
    if (i < size * size) {
        data[i] = wall_temp;
    }
}

void generate_fire_area(bool* fire_area, int idx) {
    // generate the fire area
    int i = idx / size;
    int j = idx % size;

    fire_area[i * size + j] = 0;
    int a = 0, b = 0, r2 = 0;

    float fire1_r2 = fire_size * fire_size;
    a = i - size / 2;
    b = j - size / 2;
    r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
    if (r2 < fire1_r2)
        fire_area[i * size + j] = 1;

    float fire2_r2 = (fire_size / 2) * (fire_size / 2);
    a = i - 1 * size / 3;
    b = j - 1 * size / 3;
    r2 = a * a + b * b;
    if (r2 < fire2_r2)
        fire_area[i * size + j] = 1;
}

void update(float* data, float* new_data, int i) {
    // update the temperature of each point, and store the result in `new_data`
    // to avoid data racing
    if (i < size * size) {
        if (!(i % size == 0 || i % size == size - 1 || i < size ||
              i > size * (size - 1))) {
            float up = data[i - size];
            float down = data[i + size];
            float left = data[i - 1];
            float right = data[i + 1];
            float new_val = (up + down + left + right) / 4;
            new_data[i] = new_val;
        }
    }
}

void maintain_fire(float* data, bool* fire_area, int i) {
    // maintain the temperature of fire
    if (i < size * size)
        if (fire_area[i])
            data[i] = fire_temp;
}

void maintain_wall(float* data, int i) {
    // TODO: maintain the temperature of the wall
    if (i < size * size)
        if (i % size == 0 || i % size == size - 1 || i < size ||
            i > size * (size - 1))
            data[i] = wall_temp;
}

#ifdef GUI
void data2pixels(float* data, GLubyte* pixels, int idx) {
    // convert rawdata (large, size^2) to pixels (small, resolution^2) for
    // faster rendering speed
    float factor_data_pixel = (float)size / resolution;
    float factor_temp_color = (float)255 / fire_temp;

    int x = idx / resolution;
    int y = idx % resolution;
    int idx_pixel = idx * 3;
    int x_raw = x * factor_data_pixel;
    int y_raw = y * factor_data_pixel;
    int idx_raw = y_raw * size + x_raw;
    float temp = data[idx_raw];
    int color = ((int)temp / 5 * 5) * factor_temp_color;
    pixels[idx_pixel] = color;
    pixels[idx_pixel + 1] = 255 - color;
    pixels[idx_pixel + 2] = 255 - color;
}

void plot(GLubyte* pixels) {
// visualize temprature distribution
#ifdef GUI
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(resolution, resolution, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glutSwapBuffers();
#endif
}
#endif

void master() {
    float* data_odd;
    float* data_even;
    bool* fire_area;

#ifdef GUI
    GLubyte* pixels = new GLubyte[resolution * resolution * 3];
#endif

    data_odd = new float[size * size];
    data_even = new float[size * size];

    fire_area = new bool[size * size];
    omp_set_num_threads(n_omp_threads);
#pragma omp parallel for
    for (int i = 0; i < size * size; ++i) {
        generate_fire_area(fire_area, i);
    }
    omp_set_num_threads(n_omp_threads);
#pragma omp parallel for
    for (int i = 0; i < size * size; ++i) {
        initialize(data_odd, i);
    }

    int count = 1;
    double total_time = 0;

    while (true) {
        std::chrono::high_resolution_clock::time_point t1 =
            std::chrono::high_resolution_clock::now();
        omp_set_num_threads(n_omp_threads);
#pragma omp parallel for
        for (int i = 0; i < size * size; ++i) {
            if (count % 2 == 1) {
                update(data_odd, data_even, i);
                maintain_fire(data_even, fire_area, i);
                maintain_wall(data_even, i);
            } else {
                update(data_even, data_odd, i);
                maintain_fire(data_odd, fire_area, i);
                maintain_wall(data_odd, i);
            }
        }

        std::chrono::high_resolution_clock::time_point t2 =
            std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;

#ifdef GUI
#pragma omp parallel for
        for (int i = 0; i < resolution * resolution; ++i) {
            if (count % 2 == 1) {
                data2pixels(data_even, pixels, i);
            } else {
                data2pixels(data_odd, pixels, i);
            }
        }

        plot(pixels);
#endif
        if (count >= n_iter)
            break;
    }

    printf(
        "Converge after %d iterations, elapsed time: %.6f, average computation "
        "time: %.6f\n",
        count - 1, total_time, (double)total_time / (count - 1));

    delete[] data_odd;
    delete[] data_even;
    delete[] fire_area;

#ifdef GUI
    delete[] pixels;
#endif
}

int main(int argc, char* argv[]) {
    size = atoi(argv[1]);
    n_omp_threads = atoi(argv[2]);
    n_iter = atoi(argv[3]);

#ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(resolution, resolution);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, resolution, 0, resolution);
#endif
    master();

    printf("Student ID: 120090453\n");  // replace it with your student id
    printf("Name: Haonan XUE\n");       // replace it with your name
    printf(
        "Assignment 4: Heat Distribution Simulation OpenMP Implementation\n");

    return 0;
}
