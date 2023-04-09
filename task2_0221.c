/*programing task2*/
/*2次元上の3角形2つを最も良く重ね合わせる。なお、2つの3角形は任意で異なる形状とし、プログラム実行時に入力できるようにする。*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

const double learning_rate_axis = 0.003;
const double learning_rate_angle = 0.001;

void read_file(double[][2], double[][2]);
double cross_product(double, double, double, double, double, double);
void moved_tri_init(double[][2], double[][2], double[][2], int);
void rotate_tri(double[][2], double[][2], double);
double cal_area(double[][2], double[][2]);
double differential_axis_x(double[][2], double[][2]);
double differential_axis_y(double[][2], double[][2]);
double differential_angle(double[][2], double[][2]);
void plot_result(double[][2], double[][2], double[][2]);

int main()
{
    double tri1[3][2];
    double tri2[3][2];
    double moved_tri2[3][2];
    double max_tri[3][2]; // coords of a triangle with the largest area
    double max_area = 0.0;

    read_file(tri1, tri2);
    if (cross_product(tri1[0][0], tri1[0][1], tri1[1][0], tri1[1][1], tri1[2][0], tri1[2][1]) >= 0 || cross_product(tri2[0][0], tri2[0][1], tri2[1][0], tri2[1][1], tri2[2][0], tri2[2][1]) >= 0)
    {
        printf("INPUT ERROR\nthe coordinates enterd are in the wrong order\n");
        return 0;
    }

    for (int i = 0; i < 3; i++)
    {
        moved_tri2[i][0] = tri2[i][0];
        moved_tri2[i][1] = tri2[i][1];
    }

    for (int loop_i = 0; loop_i < 9; loop_i++)
    {
        int loop_counter = 0;
        int last_change = 0;
        double area = 0.0;
        double in_loop_max_area = 0.0;
        double in_loop_max_tri[3][2];

        // initialization
        moved_tri_init(tri1, tri2, moved_tri2, loop_i);

        // serch for the local optimal solutions
        while (loop_counter - last_change < 3)
        {
            loop_counter++;

            double move_angle = learning_rate_angle * differential_angle(tri1, moved_tri2);
            rotate_tri(moved_tri2, moved_tri2, move_angle);

            double diff_x = differential_axis_x(tri1, moved_tri2);
            double diff_y = differential_axis_y(tri1, moved_tri2);

            for (int loop_j = 0; loop_j < 3; loop_j++)
            {
                moved_tri2[loop_j][0] += learning_rate_axis * diff_x;
                moved_tri2[loop_j][1] += learning_rate_axis * diff_y;
            }

            area = cal_area(tri1, moved_tri2);

            if (area > in_loop_max_area)
            {
                in_loop_max_area = area;
                last_change = loop_counter;
                for (int loop_j = 0; loop_j < 3; loop_j++)
                {
                    in_loop_max_tri[loop_j][0] = moved_tri2[loop_j][0];
                    in_loop_max_tri[loop_j][1] = moved_tri2[loop_j][1];
                }
            }
        }

        printf("in_loop_max_area[%d]: %lf\n", loop_i, in_loop_max_area);
        printf("loop_counter: %d\n", loop_counter);

        if (in_loop_max_area > max_area)
        {
            max_area = in_loop_max_area;
            for (int loop_j = 0; loop_j < 3; loop_j++)
            {
                max_tri[loop_j][0] = in_loop_max_tri[loop_j][0];
                max_tri[loop_j][1] = in_loop_max_tri[loop_j][1];
            }
        }
    }

    double sec = (double)clock() / CLOCKS_PER_SEC;

    // print result
    printf("\n-----input data-----\n");
    printf("tri1:(%lf,%lf)(%lf,%lf)(%lf,%lf)\n", tri1[0][0], tri1[0][1], tri1[1][0], tri1[1][1], tri1[2][0], tri1[2][1]);
    printf("tri2:(%lf,%lf)(%lf,%lf)(%lf,%lf)\n", tri2[0][0], tri2[0][1], tri2[1][0], tri2[1][1], tri2[2][0], tri2[2][1]);
    printf("-----output data-----\n");
    printf("tri2:(%lf,%lf)(%lf,%lf)(%lf,%lf)\n", max_tri[0][0], max_tri[0][1], max_tri[1][0], max_tri[1][1], max_tri[2][0], max_tri[2][1]);
    printf("time: %lfsec\n", sec);
    printf("max area: %lf\n", max_area);
    if (max_area != cal_area(tri1, max_tri))
    {
        printf("failed\n");
    }
    plot_result(tri1, tri2, max_tri);

    return 0;
}

void read_file(double tri1[3][2], double tri2[3][2])
{
    FILE *input;
    input = fopen("input.txt", "r");
    fscanf(input, " (%lf,%lf)(%lf,%lf)(%lf,%lf)", &tri1[0][0], &tri1[0][1], &tri1[1][0], &tri1[1][1], &tri1[2][0], &tri1[2][1]);
    fscanf(input, " (%lf,%lf)(%lf,%lf)(%lf,%lf)", &tri2[0][0], &tri2[0][1], &tri2[1][0], &tri2[1][1], &tri2[2][0], &tri2[2][1]);
    fclose(input);

    FILE *output;
    output = fopen("output1.txt", "w");
    for (int i = 0; i < 3; i++)
    {
        fprintf(output, "%lf %lf\n", tri1[i][0], tri1[i][1]);
    }
    fprintf(output, "%lf %lf\n", tri1[0][0], tri1[0][1]);
    fclose(output);
    output = fopen("output2.txt", "w");
    for (int i = 0; i < 3; i++)
    {
        fprintf(output, "%lf %lf\n", tri2[i][0], tri2[i][1]);
    }
    fprintf(output, "%lf %lf\n", tri2[0][0], tri2[0][1]);
    fclose(output);
}

void plot_result(double tri1[3][2], double tri2[3][2], double max_tri[3][2])
{
    FILE *output;
    output = fopen("output3.txt", "w");
    for (int i = 0; i < 3; i++)
    {
        fprintf(output, "%lf %lf\n", max_tri[i][0], max_tri[i][1]);
    }
    fprintf(output, "%lf %lf\n", max_tri[0][0], max_tri[0][1]);
    fclose(output);

    FILE *gp;
    gp = popen("gnuplot", "w");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set tics font 'Arial,15'\n");
    fprintf(gp, "set key font 'Arial,15'\n");
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "plot 'output2.txt' using 1:2 w l linewidth 5 lc rgb 'grey70' title 'input tri2' \n");
    fprintf(gp, "replot 'output1.txt' using 1:2 w l linewidth 5 lc rgb 'steelblue' title 'input tri1'\n");
    fprintf(gp, "replot 'output3.txt' using 1:2 w l linewidth 5 lc rgb 'red' title 'output tri2' \n");
    fprintf(gp, "set output 'output.png'\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "replot\n");
    pclose(gp);
}

void cal_intersection(double tri1[3][2], double moved_tri2[3][2], double point_set[12][2], int *point_set_num)
{
    for (int loop_i = 0; loop_i < 3; loop_i++)
    {
        // loop_c1 : 0 1 2
        // loop_c2 : 1 2 0
        int loop_c1 = loop_i;
        int loop_c2 = loop_i + 1;
        if (loop_i == 2)
        {
            loop_c2 = 0;
        }

        for (int point_loop = 0; point_loop < 3; point_loop++)
        {
            // point_c1 : 0 1 2
            // point_c2 : 1 2 0
            int point_c1 = point_loop;
            int point_c2 = point_loop + 1;
            if (point_loop == 2)
            {
                point_c2 = 0;
            }

            double tmp_n1 = (tri1[loop_c2][0] - tri1[loop_c1][0]) * (moved_tri2[point_c1][1] - tri1[loop_c1][1]) - (moved_tri2[point_c1][0] - tri1[loop_c1][0]) * (tri1[loop_c2][1] - tri1[loop_c1][1]);
            double tmp_n2 = (tri1[loop_c2][0] - tri1[loop_c1][0]) * (moved_tri2[point_c2][1] - tri1[loop_c1][1]) - (moved_tri2[point_c2][0] - tri1[loop_c1][0]) * (tri1[loop_c2][1] - tri1[loop_c1][1]);

            double tmp_n3 = (moved_tri2[point_c2][0] - moved_tri2[point_c1][0]) * (tri1[loop_c1][1] - moved_tri2[point_c1][1]) - (tri1[loop_c1][0] - moved_tri2[point_c1][0]) * (moved_tri2[point_c2][1] - moved_tri2[point_c1][1]);
            double tmp_n4 = (moved_tri2[point_c2][0] - moved_tri2[point_c1][0]) * (tri1[loop_c2][1] - moved_tri2[point_c1][1]) - (tri1[loop_c2][0] - moved_tri2[point_c1][0]) * (moved_tri2[point_c2][1] - moved_tri2[point_c1][1]);

            if (tmp_n1 * tmp_n2 < 0 && tmp_n3 * tmp_n4 < 0)
            {
                double tmp = (tri1[loop_c2][1] - tri1[loop_c1][1]) * (moved_tri2[point_c2][0] - moved_tri2[point_c1][0]) - (tri1[loop_c2][0] - tri1[loop_c1][0]) * (moved_tri2[point_c2][1] - moved_tri2[point_c1][1]);

                point_set[*point_set_num][0] = (moved_tri2[point_c1][1] * moved_tri2[point_c2][0] - moved_tri2[point_c1][0] * moved_tri2[point_c2][1]) * (tri1[loop_c2][0] - tri1[loop_c1][0]) - (tri1[loop_c1][1] * tri1[loop_c2][0] - tri1[loop_c1][0] * tri1[loop_c2][1]) * (moved_tri2[point_c2][0] - moved_tri2[point_c1][0]);
                point_set[*point_set_num][1] = (moved_tri2[point_c1][1] * moved_tri2[point_c2][0] - moved_tri2[point_c1][0] * moved_tri2[point_c2][1]) * (tri1[loop_c2][1] - tri1[loop_c1][1]) - (tri1[loop_c1][1] * tri1[loop_c2][0] - tri1[loop_c1][0] * tri1[loop_c2][1]) * (moved_tri2[point_c2][1] - moved_tri2[point_c1][1]);
                point_set[*point_set_num][0] /= tmp;
                point_set[*point_set_num][1] /= tmp;
                *point_set_num += 1;
            }
        }
    }
}

// pre_tri ->  rotate(radian) -> post_tri
void rotate_tri(double pre_tri[3][2], double post_tri[3][2], double radian)
{
    double inner_center[2]; // center of inscribed circle

    for (int i = 0; i < 3; i++)
    {
        post_tri[i][0] = pre_tri[i][0];
        post_tri[i][1] = pre_tri[i][1];
    }

    // temporary variable to calculate inner_center
    double tmp_s = sqrt(pow(post_tri[2][0] - post_tri[1][0], 2.0) + pow(post_tri[2][1] - post_tri[1][1], 2.0));
    double tmp_t = sqrt(pow(post_tri[0][0] - post_tri[2][0], 2.0) + pow(post_tri[0][1] - post_tri[2][1], 2.0));
    double tmp_u = sqrt(pow(post_tri[1][0] - post_tri[0][0], 2.0) + pow(post_tri[1][1] - post_tri[0][1], 2.0));

    inner_center[0] = (tmp_s * post_tri[0][0] + tmp_t * post_tri[1][0] + tmp_u * post_tri[2][0]) / (tmp_s + tmp_t + tmp_u);
    inner_center[1] = (tmp_s * post_tri[0][1] + tmp_t * post_tri[1][1] + tmp_u * post_tri[2][1]) / (tmp_s + tmp_t + tmp_u);

    for (int i = 0; i < 3; i++)
    {
        post_tri[i][0] -= inner_center[0];
        post_tri[i][1] -= inner_center[1];
    }

    for (int sc = 0; sc < 3; sc++)
    {
        double tmp_x = post_tri[sc][0] * cos(radian) - post_tri[sc][1] * sin(radian);
        double tmp_y = post_tri[sc][0] * sin(radian) + post_tri[sc][1] * cos(radian);
        post_tri[sc][0] = tmp_x;
        post_tri[sc][1] = tmp_y;
    }

    for (int i = 0; i < 3; i++)
    {
        post_tri[i][0] += inner_center[0];
        post_tri[i][1] += inner_center[1];
    }
}

double cross_product(double point1_x, double point1_y, double point2_x, double point2_y, double point3_x, double point3_y)
{
    double vec[2][2];
    vec[0][0] = point2_x - point1_x;
    vec[0][1] = point2_y - point1_y;
    vec[1][0] = point3_x - point2_x;
    vec[1][1] = point3_y - point2_y;
    double cross_p = vec[0][0] * vec[1][1] - vec[0][1] * vec[1][0];
    return cross_p;
}

bool judge_io(double judge_point_x, double judge_point_y, double tri[][2])
{
    double cp1 = cross_product(tri[2][0], tri[2][1], tri[0][0], tri[0][1], judge_point_x, judge_point_y);
    double cp2 = cross_product(tri[0][0], tri[0][1], tri[1][0], tri[1][1], judge_point_x, judge_point_y);
    double cp3 = cross_product(tri[1][0], tri[1][1], tri[2][0], tri[2][1], judge_point_x, judge_point_y);
    if ((cp1 >= 0 && cp2 >= 0 && cp3 >= 0) || (cp1 <= 0 && cp2 <= 0 && cp3 <= 0))
    {
        return true; // inside triangle
    }
    else
    {
        return false; // outside triangle
    }
}

bool judge_new(double judge_point_x, double judge_point_y, double point_set[][2], int point_set_num)
{
    for (int loop_i = 0; loop_i < point_set_num; loop_i++)
    {
        if (judge_point_x == point_set[loop_i][0] && judge_point_y == point_set[loop_i][1])
        {
            return false; // overlapping
        }
    }
    return true; // new point
}

void add_vertex(double point_set[12][2], int *point_set_num, double tri1[][2], double moved_tri2[][2])
{
    // If each point is inside a triangle, add it to the array
    for (int i = 0; i < 3; i++)
    {
        if (judge_io(tri1[i][0], tri1[i][1], moved_tri2) == true && judge_new(tri1[i][0], tri1[i][1], point_set, *point_set_num))
        {
            point_set[*point_set_num][0] = tri1[i][0];
            point_set[*point_set_num][1] = tri1[i][1];
            *point_set_num += 1;
        }
        if (judge_io(moved_tri2[i][0], moved_tri2[i][1], tri1) == true && judge_new(moved_tri2[i][0], moved_tri2[i][1], point_set, *point_set_num))
        {
            point_set[*point_set_num][0] = moved_tri2[i][0];
            point_set[*point_set_num][1] = moved_tri2[i][1];
            *point_set_num += 1;
        }
    }
}

void argument_sort(double point_set[][2], int point_set_num)
{
    double min_y[2]; // point with the smallest y
    min_y[0] = point_set[0][0];
    min_y[1] = point_set[0][1];
    for (int i = 1; i < point_set_num; i++)
    {
        if (point_set[i][1] < min_y[1])
        {
            min_y[0] = point_set[i][0];
            min_y[1] = point_set[i][1];
        }
    }

    for (int i = 0; i < point_set_num; i++)
    {
        point_set[i][0] -= min_y[0];
        point_set[i][1] -= min_y[1];
    }
    for (int sc1 = 0; sc1 < point_set_num; sc1++)
    {
        for (int sc2 = sc1 + 1; sc2 < point_set_num; sc2++)
        {
            if (atan2(point_set[sc1][1], point_set[sc1][0]) > atan2(point_set[sc2][1], point_set[sc2][0]))
            {
                double tmp1 = point_set[sc1][0];
                double tmp2 = point_set[sc1][1];
                point_set[sc1][0] = point_set[sc2][0];
                point_set[sc1][1] = point_set[sc2][1];
                point_set[sc2][0] = tmp1;
                point_set[sc2][1] = tmp2;
            }
        }
    }
    for (int i = 0; i < point_set_num; i++)
    {
        point_set[i][0] += min_y[0];
        point_set[i][1] += min_y[1];
    }
}

double coordinates_method(double polygon[][2], int in_num)
{
    double area = 0.0;
    for (int ac = 0; ac < in_num - 1; ac++)
    {
        area += polygon[ac][0] * polygon[ac + 1][1] - polygon[ac + 1][0] * polygon[ac][1];
    }
    area += polygon[in_num - 1][0] * polygon[0][1] - polygon[0][0] * polygon[in_num - 1][1];
    area *= 0.5;
    area = fabs(area);
    return area;
}

double cal_area(double tri1[][2], double moved_tri2[][2])
{
    double point_set[12][2]; // coordinates composing a polygon
    int point_set_num = 0;

    cal_intersection(tri1, moved_tri2, point_set, &point_set_num);

    add_vertex(point_set, &point_set_num, tri1, moved_tri2);

    argument_sort(point_set, point_set_num);

    double area = coordinates_method(point_set, point_set_num);

    return area;
}

double differential_axis_x(double tri1[][2], double moved_tri2[][2])
{
    const double difference_value_axis = 0.01;
    double moved_tri2_minus[3][2];
    double moved_tri2_plus[3][2];

    for (int dc = 0; dc < 3; dc++)
    {
        moved_tri2_minus[dc][0] = moved_tri2[dc][0] - difference_value_axis;
        moved_tri2_minus[dc][1] = moved_tri2[dc][1];
        moved_tri2_plus[dc][0] = moved_tri2[dc][0] + difference_value_axis;
        moved_tri2_plus[dc][1] = moved_tri2[dc][1];
    }

    double minus_area = cal_area(tri1, moved_tri2_minus);
    double plus_area = cal_area(tri1, moved_tri2_plus);

    double differential = (plus_area - minus_area) / (difference_value_axis * 2);

    return differential;
}

double differential_axis_y(double tri1[][2], double moved_tri2[][2])
{
    const double difference_value_axis = 0.01;
    double moved_tri2_minus[3][2];
    double moved_tri2_plus[3][2];

    for (int dc = 0; dc < 3; dc++)
    {
        moved_tri2_minus[dc][0] = moved_tri2[dc][0];
        moved_tri2_minus[dc][1] = moved_tri2[dc][1] - difference_value_axis;
        moved_tri2_plus[dc][0] = moved_tri2[dc][0];
        moved_tri2_plus[dc][1] = moved_tri2[dc][1] + difference_value_axis;
    }

    double minus_area = cal_area(tri1, moved_tri2_minus);
    double plus_area = cal_area(tri1, moved_tri2_plus);

    double differential = (plus_area - minus_area) / (difference_value_axis * 2);

    return differential;
}

double differential_angle(double tri1[][2], double moved_tri2[][2])
{
    const double difference_value_angle = 0.01;
    double moved_tri2_minus[3][2];
    double moved_tri2_plus[3][2];

    rotate_tri(moved_tri2, moved_tri2_minus, -difference_value_angle);
    rotate_tri(moved_tri2, moved_tri2_plus, difference_value_angle);

    double minus_area = cal_area(tri1, moved_tri2_minus);
    double plus_area = cal_area(tri1, moved_tri2_plus);

    double differential = (plus_area - minus_area) / (difference_value_angle * 2);

    return differential;
}

void moved_tri_init(double tri1[][2], double tri2[][2], double moved_tri2[][2], int loop_num)
{
    // initialization moved_tri2

    double vec_tri1[2];
    double vec_moved_tri2[2];
    double tri1_center_point[2];
    double moved_tri2_center_point[2];

    // loop : 0  1  2  3  4  5  6  7  8
    // num1 : 0  0  0  1  1  1  2  2  2
    // num2 : 0  1  2  0  1  2  0  1  2

    int num1 = loop_num / 3;
    int num2 = loop_num % 3;

    for (int sc = 0; sc < 3; sc++)
    {
        moved_tri2[sc][0] = tri2[sc][0];
        moved_tri2[sc][1] = tri2[sc][1];
    }

    switch (num1)
    {
    case 0:
        vec_tri1[0] = tri1[1][0] - tri1[0][0];
        vec_tri1[1] = tri1[1][1] - tri1[0][1];
        tri1_center_point[0] = (tri1[0][0] + tri1[1][0]) / 2;
        tri1_center_point[1] = (tri1[0][1] + tri1[1][1]) / 2;
        break;
    case 1:
        vec_tri1[0] = tri1[2][0] - tri1[1][0];
        vec_tri1[1] = tri1[2][1] - tri1[1][1];
        tri1_center_point[0] = (tri1[1][0] + tri1[2][0]) / 2;
        tri1_center_point[1] = (tri1[1][1] + tri1[2][1]) / 2;
        break;
    case 2:
        vec_tri1[0] = tri1[0][0] - tri1[2][0];
        vec_tri1[1] = tri1[0][1] - tri1[2][1];
        tri1_center_point[0] = (tri1[2][0] + tri1[0][0]) / 2;
        tri1_center_point[1] = (tri1[2][1] + tri1[0][1]) / 2;
        break;
    }

    switch (num2)
    {
    case 0:
        vec_moved_tri2[0] = moved_tri2[1][0] - moved_tri2[0][0];
        vec_moved_tri2[1] = moved_tri2[1][1] - moved_tri2[0][1];
        break;
    case 1:
        vec_moved_tri2[0] = moved_tri2[2][0] - moved_tri2[1][0];
        vec_moved_tri2[1] = moved_tri2[2][1] - moved_tri2[1][1];
        break;
    case 2:
        vec_moved_tri2[0] = moved_tri2[0][0] - moved_tri2[2][0];
        vec_moved_tri2[1] = moved_tri2[0][1] - moved_tri2[2][1];
        break;
    }

    double move_value_angle = atan2(vec_tri1[1], vec_tri1[0]) - atan2(vec_moved_tri2[1], vec_moved_tri2[0]);

    for (int i = 0; i < 3; i++)
    {
        moved_tri2[i][0] -= tri2[num2][0];
        moved_tri2[i][1] -= tri2[num2][1];
    }

    for (int sc = 0; sc < 3; sc++)
    {
        double rx = moved_tri2[sc][0] * cos(move_value_angle) - moved_tri2[sc][1] * sin(move_value_angle);
        double ry = moved_tri2[sc][0] * sin(move_value_angle) + moved_tri2[sc][1] * cos(move_value_angle);
        moved_tri2[sc][0] = rx;
        moved_tri2[sc][1] = ry;
    }

    switch (num2)
    {
    case 0:
        moved_tri2_center_point[0] = (moved_tri2[0][0] + moved_tri2[1][0]) / 2;
        moved_tri2_center_point[1] = (moved_tri2[0][1] + moved_tri2[1][1]) / 2;
        break;
    case 1:
        moved_tri2_center_point[0] = (moved_tri2[1][0] + moved_tri2[2][0]) / 2;
        moved_tri2_center_point[1] = (moved_tri2[1][1] + moved_tri2[2][1]) / 2;
        break;
    case 2:
        moved_tri2_center_point[0] = (moved_tri2[2][0] + moved_tri2[0][0]) / 2;
        moved_tri2_center_point[1] = (moved_tri2[2][1] + moved_tri2[0][1]) / 2;
        break;
    }

    double move_value_x = tri1_center_point[0] - moved_tri2_center_point[0];
    double move_value_y = tri1_center_point[1] - moved_tri2_center_point[1];

    for (int sc = 0; sc < 3; sc++)
    {
        moved_tri2[sc][0] += (move_value_x);
        moved_tri2[sc][1] += (move_value_y);
    }
}
