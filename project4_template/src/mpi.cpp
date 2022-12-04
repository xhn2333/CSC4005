#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"

int size;  // problem size

int my_rank;
int world_size;

void initialize(float* data) {
    // intialize the temperature distribution
    int len = size * size;
    for (int i = 0; i < len; i++) {
        data[i] = wall_temp;
    }
}

void generate_fire_area(bool* fire_area) {
    // generate the fire area
    int len = size * size;
    for (int i = 0; i < len; i++) {
        fire_area[i] = 0;
    }

    float fire1_r2 = fire_size * fire_size;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int a = i - size / 2;
            int b = j - size / 2;
            int r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
            if (r2 < fire1_r2)
                fire_area[i * size + j] = 1;
        }
    }

    float fire2_r2 = (fire_size / 2) * (fire_size / 2);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int a = i - 1 * size / 3;
            int b = j - 1 * size / 3;
            int r2 = a * a + b * b;
            if (r2 < fire2_r2)
                fire_area[i * size + j] = 1;
        }
    }
}

void update(float* data, float* new_data, int begin, int end) {
    // TODO: update the temperature of each point, and store the result in
    // `new_data` to avoid data racing

    for (int i = begin; i < end; ++i)
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

void maintain_fire(float* data, bool* fire_area, int begin, int end) {
    // TODO: maintain the temperature of fire
    for (int i = begin; i < end; i++) {
        if (fire_area[i])
            data[i] = fire_temp;
    }
}

void maintain_wall(float* data, int begin, int end) {
    // TODO: maintain the temperature of the wall
    for (int i = begin; i < end; i++) {
        if (i % size == 0 || i % size == size - 1 || i < size ||
            i > size * (size - 1))
            data[i] = wall_temp;
    }
}

#ifdef GUI
void data2pixels(float* data, GLubyte* pixels, int begin, int end) {
    // convert rawdata (large, size^2) to pixels (small, resolution^2) for
    // faster rendering speed
    float factor_data_pixel = (float)size / resolution;
    float factor_temp_color = (float)255 / fire_temp;
    for (int x = 0; x < resolution; x++) {
        for (int y = 0; y < resolution; y++) {
            int idx = x * resolution + y;
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
    }
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

int main(int argc, char* argv[]) {
    size = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    float* data_odd = new float[size * size];
    float* data_even = new float[size * size];
    bool* fire_area = new bool[size * size];

    if (my_rank == 0) {
#ifdef GUI
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
        glutInitWindowPosition(0, 0);
        glutInitWindowSize(resolution, resolution);
        glutCreateWindow(
            "Heat Distribution Simulation Sequential Implementation");
        gluOrtho2D(0, resolution, 0, resolution);
#endif
        initialize(data_odd);
        generate_fire_area(fire_area);
    }

    int my_begin = size * size * my_rank / world_size;
    int my_end = size * size * (my_rank + 1) / world_size;

#ifdef GUI
    GLubyte* pixels;
    pixels = new GLubyte[resolution * resolution * 3];
#endif

    int count = 1;
    double total_time = 0;

    int* sendcounts = new int[world_size];
    int* displacements = new int[world_size];
    int n_element = (size * size - 1) / world_size + 1;
    for (int i = 0; i < world_size; ++i) {
        displacements[i] = i * n_element;
        sendcounts[i] = n_element < size * size - displacements[i]
                            ? n_element
                            : size * size - displacements[i];
    }

    // TODO: Send initial distribution to each slave process

    MPI_Bcast(data_even, size * size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data_odd, size * size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(fire_area, size * size, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    while (true) {
        std::chrono::high_resolution_clock::time_point t1 =
            std::chrono::high_resolution_clock::now();

        // TODO: Computation of my part
        update(data_odd, data_even, my_begin, my_end);
        maintain_fire(data_even, fire_area, my_begin, my_end);
        maintain_wall(data_even, my_begin, my_end);
        MPI_Allgatherv(data_even + displacements[my_rank], sendcounts[my_rank],
                       MPI_FLOAT, data_odd, sendcounts, displacements,
                       MPI_FLOAT, MPI_COMM_WORLD);

        // TODO: Send border row to neighbours

        std::chrono::high_resolution_clock::time_point t2 =
            std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;

#ifdef GUI
        if (my_rank == 0) {
            // TODO: Gather pixels of slave processes
            data2pixels(data_odd, pixels, my_begin, my_end);
            plot(pixels);
        }
#endif
    }

    delete[] data_odd;
    delete[] data_even;
    delete[] fire_area;

#ifdef GUI
    delete[] pixels;
#endif

    if (my_rank == 0) {
        printf("Student ID: 120090453\n");  // replace it with your student id
        printf("Name: Haonan XUE\n");       // replace it with your name
        printf(
            "Assignment 4: Heat Distribution Simulation MPI Implementation\n");
    }

    MPI_Finalize();

    return 0;
}
