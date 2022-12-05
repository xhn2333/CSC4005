#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"

int block_size = 512;     // cuda thread block size
__device__ int size_gpu;  // problem size
int size_cpu;
int n_iter;

__global__ void initialize(float* data) {
    // TODO: intialize the temperature distribution (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < size_gpu * size_gpu) {
        data[i] = wall_temp;
    }
}

__global__ void generate_fire_area(bool* fire_area) {
    // TODO: generate the fire area (in parallelized way)
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    int i = idx / size_gpu;
    int j = idx % size_gpu;

    fire_area[i * size_gpu + j] = 0;
    int a = 0, b = 0, r2 = 0;

    float fire1_r2 = fire_size * fire_size;
    a = i - size_gpu / 2;
    b = j - size_gpu / 2;
    r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
    if (r2 < fire1_r2)
        fire_area[i * size_gpu + j] = 1;

    float fire2_r2 = (fire_size / 2) * (fire_size / 2);
    a = i - 1 * size_gpu / 3;
    b = j - 1 * size_gpu / 3;
    r2 = a * a + b * b;
    if (r2 < fire2_r2)
        fire_area[i * size_gpu + j] = 1;
}

__global__ void update(float* data, float* new_data) {
    // TODO: update temperature for each point  (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < size_gpu * size_gpu) {
        if (!(i % size_gpu == 0 || i % size_gpu == size_gpu - 1 ||
              i < size_gpu || i > size_gpu * (size_gpu - 1))) {
            float up = data[i - size_gpu];
            float down = data[i + size_gpu];
            float left = data[i - 1];
            float right = data[i + 1];
            float new_val = (up + down + left + right) / 4;
            new_data[i] = new_val;
        }
    }
}

__global__ void maintain_wall(float* data) {
    // TODO: maintain the temperature of the wall (sequential is enough)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < size_gpu * size_gpu) {
        if (i % size_gpu == 0 || i % size_gpu == size_gpu - 1 || i < size_gpu ||
            i > size_gpu * (size_gpu - 1))
            data[i] = wall_temp;
    }
}

__global__ void maintain_fire(float* data, bool* fire_area) {
    // TODO: maintain the temperature of the fire (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < size_gpu * size_gpu) {
        if (fire_area[i])
            data[i] = fire_temp;
    }
}

#ifdef GUI
__global__ void data2pixels(float* data, GLubyte* pixels) {
    // TODO: convert rawdata (large, size^2) to pixels (small, resolution^2) for
    // faster rendering speed (in parallelized way)
    float factor_data_pixel = (float)size_gpu / resolution;
    float factor_temp_color = (float)255 / fire_temp;

    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    int x = idx / resolution;
    int y = idx % resolution;
    int idx_pixel = idx * 3;
    int x_raw = x * factor_data_pixel;
    int y_raw = y * factor_data_pixel;
    int idx_raw = y_raw * size_gpu + x_raw;
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

__global__ void warmup() {}

void master() {
    float* data_odd;
    float* data_even;
    bool* fire_area;

    clock_t clock_start;
    clock_t clock_end;

    cudaMalloc(&data_odd, size_cpu * size_cpu * sizeof(float));
    cudaMalloc(&data_even, size_cpu * size_cpu * sizeof(float));
    cudaMalloc(&fire_area, size_cpu * size_cpu * sizeof(bool));

#ifdef GUI
    GLubyte* pixels;
    GLubyte* host_pixels;
    host_pixels = new GLubyte[resolution * resolution * 3];
    cudaMalloc(&pixels, resolution * resolution * 3 * sizeof(GLubyte));
#endif

    int n_block_size = size_cpu * size_cpu / block_size + 1;
    int n_block_resolution = resolution * resolution / block_size + 1;

    initialize<<<n_block_size, block_size>>>(data_odd);
    generate_fire_area<<<n_block_size, block_size>>>(fire_area);

    int count = 1;
    double total_time = 0;

    while (true) {
        // std::chrono::high_resolution_clock::time_point t1 =
        // std::chrono::high_resolution_clock::now();
        clock_start = clock();
        // TODO: modify the following lines to fit your need.
        if (count % 2 == 1) {
            update<<<n_block_size, block_size>>>(data_odd, data_even);
            maintain_fire<<<n_block_size, block_size>>>(data_even, fire_area);
            maintain_wall<<<n_block_size, block_size>>>(data_even);
        } else {
            update<<<n_block_size, block_size>>>(data_even, data_odd);
            maintain_fire<<<n_block_size, block_size>>>(data_odd, fire_area);
            maintain_wall<<<n_block_size, block_size>>>(data_odd);
        }
        clock_end = clock();
        // std::chrono::high_resolution_clock::time_point t2 =
        // std::chrono::high_resolution_clock::now();
        double this_time =
            double(clock_end - clock_start) / double(CLOCKS_PER_SEC);
        // std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;

#ifdef GUI
        if (count % 2 == 1) {
            data2pixels<<<n_block_resolution, block_size>>>(data_even, pixels);
        } else {
            data2pixels<<<n_block_resolution, block_size>>>(data_odd, pixels);
        }
        cudaMemcpy(host_pixels, pixels,
                   resolution * resolution * 3 * sizeof(GLubyte),
                   cudaMemcpyDeviceToHost);
        plot(host_pixels);
#endif
        if (count >= n_iter)
            break;
    }

    printf(
        "Converge after %d iterations, elapsed time: %.6f, average computation "
        "time: %.6f\n",
        count - 1, total_time, (double)total_time / (count - 1));

    cudaFree(data_odd);
    cudaFree(data_even);
    cudaFree(fire_area);

#ifdef GUI
    cudaFree(pixels);
    delete[] host_pixels;
#endif
}

int main(int argc, char* argv[]) {
    size_cpu = atoi(argv[1]);
    block_size = atoi(argv[2]);
    n_iter = atoi(argv[3]);
    cudaMemcpyToSymbol(size_gpu, &size_cpu, sizeof(int));

#ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(resolution, resolution);
    glutCreateWindow("Heat Distribution Simulation Sequential Implementation");
    gluOrtho2D(0, resolution, 0, resolution);
#endif

    master();

    printf("Student ID: 120090453\n");  // replace it with your student id
    printf("Name: Haonan XUE\n");       // replace it with your name
    printf("Assignment 4: Heat Distribution CUDA Implementation\n");

    return 0;
}
