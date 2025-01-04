#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double* generate_uniform(int n, double a, double b);
double* generate_exponential(int n, double lambda);
double* generate_normal(int n, double mu, double sigma);
void calculate_moments(double* samples, int n, const char* distribution_name, double theoretical_mean, double theoretical_variance, double theoretical_skewness, double theoretical_kurtosis);
void export_to_csv(const char* filename, double* samples, int n);

int main() {
    srand(time(NULL));

    int n = 10000;

    double* uniform_samples = generate_uniform(n, 0.0, 1.0);
    double* exponential_samples = generate_exponential(n, 1.0);
    double* normal_samples = generate_normal(n, 0.0, 1.0);

    calculate_moments(uniform_samples, n, "Uniform", 0.5, 1.0 / 12.0, 0.0, -1.2);
    calculate_moments(exponential_samples, n, "Exponential", 1.0, 1.0, 2.0, 6.0);
    calculate_moments(normal_samples, n, "Normal", 0.0, 1.0, 0.0, 0.0);

    //  to CSV 
    export_to_csv("uniform_samples.csv", uniform_samples, n);
    export_to_csv("exponential_samples.csv", exponential_samples, n);
    export_to_csv("normal_samples.csv", normal_samples, n);

    printf("\nData exported to CSV files: uniform_samples.csv, exponential_samples.csv, normal_samples.csv\n");

    // Free allocated memory
    free(uniform_samples);
    free(exponential_samples);
    free(normal_samples);

    return 0;
}

//inverse transform sampling
double* generate_uniform(int n, double a, double b) {
    double* samples = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        double u = (double)rand() / RAND_MAX;
        samples[i] = a + (b - a) * u;
    }
    return samples;
}

//exponential distribution 
double* generate_exponential(int n, double lambda) {
    double* samples = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        double u = (double)rand() / RAND_MAX;
        samples[i] = -log(1 - u) / lambda;
    }
    return samples;
}

//normal distribution using the Box-Muller method
double* generate_normal(int n, double mu, double sigma) {
    double* samples = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i += 2) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
        samples[i] = mu + sigma * z0;
        if (i + 1 < n) samples[i + 1] = mu + sigma * z1;
    }
    return samples;
}

// Calculate moments 
void calculate_moments(double* samples, int n, const char* distribution_name, double theoretical_mean, double theoretical_variance, double theoretical_skewness, double theoretical_kurtosis) {
    double sum = 0.0, sum_squared = 0.0, sum_cubed = 0.0, sum_quartic = 0.0;
    for (int i = 0; i < n; i++) {
        sum += samples[i];
        sum_squared += samples[i] * samples[i];
        sum_cubed += samples[i] * samples[i] * samples[i];
        sum_quartic += samples[i] * samples[i] * samples[i] * samples[i];
    }

    double mean = sum / n;
    double variance = (sum_squared / n) - (mean * mean);
    double skewness = ((sum_cubed / n) - 3 * mean * variance - mean * mean * mean) / pow(variance, 1.5);
    double kurtosis = ((sum_quartic / n) - 4 * mean * (sum_cubed / n) + 6 * mean * mean * variance + 3 * mean * mean * mean * mean) / (variance * variance) - 3;

    printf("\n%s Distribution:\n", distribution_name);
    printf("Sample Mean: %.4f, Theoretical Mean: %.4f\n", mean, theoretical_mean);
    printf("Sample Variance: %.4f, Theoretical Variance: %.4f\n", variance, theoretical_variance);
    printf("Sample Skewness: %.4f, Theoretical Skewness: %.4f\n", skewness, theoretical_skewness);
    printf("Sample Kurtosis: %.4f, Theoretical Kurtosis: %.4f\n", kurtosis, theoretical_kurtosis);
}

// to CSV file
void export_to_csv(const char* filename, double* samples, int n) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) {
        fprintf(file, "%.6f\n", samples[i]);
    }

    fclose(file);
}
