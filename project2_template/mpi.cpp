#include <mpi.h>
#include <stdio.h>
#include <ostream>
#include "asg2.h"

int rank;
int world_size;
int num_my_element;
int *point_x, *point_y;
float* point_c;
int *global_x, *global_y;
float* global_c;

void master() {
    // TODO: procedure run in master process
    //  MPI_Scatter...
    //  MPI_Gather...
    //  the following code is not a necessary, please replace it with MPI
    //  implementation.

    // TODO END
}

void slave() {
    // TODO: procedure run in slave process
    //  MPI_Scatter...
    //  MPI_Gather...

    // TODO END
}

int main(int argc, char* argv[]) {
    if (argc == 4) {
        X_RESN = atoi(argv[1]);
        Y_RESN = atoi(argv[2]);
        max_iteration = atoi(argv[3]);
        total_size = X_RESN * Y_RESN;
    } else {
        X_RESN = 1000;
        Y_RESN = 1000;
        max_iteration = 100;
        total_size = X_RESN * Y_RESN;
    }

    /* computation part begin */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (rank == 0) {
#ifdef GUI
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
        glutInitWindowSize(500, 500);
        glutInitWindowPosition(0, 0);
        glutCreateWindow("MPI");
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glMatrixMode(GL_PROJECTION);
        gluOrtho2D(0, X_RESN, 0, Y_RESN);
        glutDisplayFunc(plot);
#endif
    }

    num_my_element = total_size / world_size;
    point_x = new int[num_my_element];
    point_y = new int[num_my_element];
    point_c = new float[num_my_element];
    global_x = new int[total_size];
    global_y = new int[total_size];
    global_c = new float[total_size];

    if (rank == 0) {
        t1 = std::chrono::high_resolution_clock::now();

        initData();
        for (int i = 0; i < total_size; ++i) {
            global_x[i] = data[i].x;
            global_y[i] = data[i].y;
            global_c[i] = data[i].color;
        }

        t2 = std::chrono::high_resolution_clock::now();
        time_span = t2 - t1;

        printf("Student ID: 120090453\n");  // replace it with your student id
        printf("Name: Haonan Xue\n");       // replace it with your name
        printf("Assignment 2 MPI\n");
        printf("Run Time: %f seconds\n", time_span.count());
        printf("Problem Size: %d * %d, %d\n", X_RESN, Y_RESN, max_iteration);
        printf("Process Number: %d\n", world_size);
    }

    MPI_Scatter(global_x, num_my_element, MPI_INT, point_x, num_my_element,
                MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(global_y, num_my_element, MPI_INT, point_y, num_my_element,
                MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(global_c, num_my_element, MPI_FLOAT, point_c, num_my_element,
                MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < num_my_element; ++i) {
        Point t;
        t.x = point_x[i];
        t.y = point_y[i];
        t.color = point_c[i];
        compute(&t);
        point_c[i] = t.color;
    }
    MPI_Gather(point_c, num_my_element, MPI_INT, global_c, num_my_element,
               MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < total_size; ++i) {
            data[i].color = global_c[i];
        }
    }
    MPI_Finalize();
    /* computation part end */

    if (rank == 0) {
#ifdef GUI
        glutMainLoop();
#endif
    }

    return 0;
}
