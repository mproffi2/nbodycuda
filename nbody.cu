#include <iostream> 
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <chrono>
#include <cuda_runtime.h>

double G = 6.674*std::pow(10,-11);
//double G = 1;

struct simulation {
  size_t nbpart;
  
  std::vector<double> mass;

  //position
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  //velocity
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;

  //force
  std::vector<double> fx;
  std::vector<double> fy;
  std::vector<double> fz;

  
  simulation(size_t nb)
    :nbpart(nb), mass(nb),
     x(nb), y(nb), z(nb),
     vx(nb), vy(nb), vz(nb),
     fx(nb), fy(nb), fz(nb) 
  {}
};


void random_init(simulation& s) {
  std::random_device rd;  
  std::mt19937 gen(rd());
  std::uniform_real_distribution dismass(0.9, 1.);
  std::normal_distribution dispos(0., 1.);
  std::normal_distribution disvel(0., 1.);

  for (size_t i = 0; i<s.nbpart; ++i) {
    s.mass[i] = dismass(gen);

    s.x[i] = dispos(gen);
    s.y[i] = dispos(gen);
    s.z[i] = dispos(gen);
    s.z[i] = 0.;
    
    s.vx[i] = disvel(gen);
    s.vy[i] = disvel(gen);
    s.vz[i] = disvel(gen);
    s.vz[i] = 0.;
    s.vx[i] = s.y[i]*1.5;
    s.vy[i] = -s.x[i]*1.5;
  }

  return;
  //normalize velocity (using normalization found on some physicis blog)
  double meanmass = 0;
  double meanmassvx = 0;
  double meanmassvy = 0;
  double meanmassvz = 0;
  for (size_t i = 0; i<s.nbpart; ++i) {
    meanmass += s.mass[i];
    meanmassvx += s.mass[i] * s.vx[i];
    meanmassvy += s.mass[i] * s.vy[i];
    meanmassvz += s.mass[i] * s.vz[i];
  }
  for (size_t i = 0; i<s.nbpart; ++i) {
    s.vx[i] -= meanmassvx/meanmass;
    s.vy[i] -= meanmassvy/meanmass;
    s.vz[i] -= meanmassvz/meanmass;
  }
  
}

void init_solar(simulation& s) {
  enum Planets {SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, MOON};
  s = simulation(10);

  // Masses in kg
  s.mass[SUN] = 1.9891 * std::pow(10, 30);
  s.mass[MERCURY] = 3.285 * std::pow(10, 23);
  s.mass[VENUS] = 4.867 * std::pow(10, 24);
  s.mass[EARTH] = 5.972 * std::pow(10, 24);
  s.mass[MARS] = 6.39 * std::pow(10, 23);
  s.mass[JUPITER] = 1.898 * std::pow(10, 27);
  s.mass[SATURN] = 5.683 * std::pow(10, 26);
  s.mass[URANUS] = 8.681 * std::pow(10, 25);
  s.mass[NEPTUNE] = 1.024 * std::pow(10, 26);
  s.mass[MOON] = 7.342 * std::pow(10, 22);

  // Positions (in meters) and velocities (in m/s)
  double AU = 1.496 * std::pow(10, 11); // Astronomical Unit

  s.x = {0, 0.39*AU, 0.72*AU, 1.0*AU, 1.52*AU, 5.20*AU, 9.58*AU, 19.22*AU, 30.05*AU, 1.0*AU + 3.844*std::pow(10, 8)};
  s.y = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  s.z = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  s.vx = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  s.vy = {0, 47870, 35020, 29780, 24130, 13070, 9680, 6800, 5430, 29780 + 1022};
  s.vz = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}

//meant to update the force that from applies on to
void update_force(simulation& s, size_t from, size_t to) {
  double softening = .1;
  double dist_sq = std::pow(s.x[from]-s.x[to],2)
    + std::pow(s.y[from]-s.y[to],2)
    + std::pow(s.z[from]-s.z[to],2);
  double F = G * s.mass[from]*s.mass[to]/(dist_sq+softening); //that the strength of the force

  //direction
  double dx = s.x[from]-s.x[to];
  double dy = s.y[from]-s.y[to];
  double dz = s.z[from]-s.z[to];
  double norm = std::sqrt(dx*dx+dy*dy+dz*dz);
  
  dx = dx/norm;
  dy = dy/norm;
  dz = dz/norm;

  //apply force
  s.fx[to] += dx*F;
  s.fy[to] += dy*F;
  s.fz[to] += dz*F;
}

void reset_force(simulation& s) {
  for (size_t i=0; i<s.nbpart; ++i) {
    s.fx[i] = 0.;
    s.fy[i] = 0.;
    s.fz[i] = 0.;
  }
}

void apply_force(simulation& s, size_t i, double dt) {
  s.vx[i] += s.fx[i]/s.mass[i]*dt;
  s.vy[i] += s.fy[i]/s.mass[i]*dt;
  s.vz[i] += s.fz[i]/s.mass[i]*dt;
}

void update_position(simulation& s, size_t i, double dt) {
  s.x[i] += s.vx[i]*dt;
  s.y[i] += s.vy[i]*dt;
  s.z[i] += s.vz[i]*dt;
}

void dump_state(simulation& s) {
  std::cout<<s.nbpart<<'\t';
  for (size_t i=0; i<s.nbpart; ++i) {
    std::cout<<s.mass[i]<<'\t';
    std::cout<<s.x[i]<<'\t'<<s.y[i]<<'\t'<<s.z[i]<<'\t';
    std::cout<<s.vx[i]<<'\t'<<s.vy[i]<<'\t'<<s.vz[i]<<'\t';
    std::cout<<s.fx[i]<<'\t'<<s.fy[i]<<'\t'<<s.fz[i]<<'\t';
  }
  std::cout<<'\n';
}

void load_from_file(simulation& s, std::string filename) {
  std::ifstream in (filename);
  size_t nbpart;
  in>>nbpart;
  s = simulation(nbpart);
  for (size_t i=0; i<s.nbpart; ++i) {
    in>>s.mass[i];
    in >>  s.x[i] >>  s.y[i] >>  s.z[i];
    in >> s.vx[i] >> s.vy[i] >> s.vz[i];
    in >> s.fx[i] >> s.fy[i] >> s.fz[i];
  }
  if (!in.good())
    throw "kaboom";
}

// softening for GPU kernel & CPU use
static const double SOFTENING = 1e-9;

// CPU-only helper that times the main loop
double run_cpu(simulation &s, double dt, size_t nbstep, size_t printevery) {
  auto tstart = std::chrono::high_resolution_clock::now();

  for (size_t step = 0; step< nbstep; step++) {
    if (step %printevery == 0)
      dump_state(s);
  
    reset_force(s);
    for (size_t i=0; i<s.nbpart; ++i)
      for (size_t j=0; j<s.nbpart; ++j)
    if (i != j)
      update_force(s, i, j);

    for (size_t i=0; i<s.nbpart; ++i) {
      apply_force(s, i, dt);
      update_position(s, i, dt);
    }
  }

  auto tend = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = tend - tstart;
  return diff.count();
}

// simple CUDA error-check helper
#define CUDA_CHECK(call) do { \
  cudaError_t err = (call); \
  if (err != cudaSuccess) { \
    fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
    exit(EXIT_FAILURE); \
  } \
} while(0)

// CUDA kernels 
// compute forces: thread i computes net force on particle i
__global__ void compute_forces_gpu(size_t n,
                                   const double* mass,
                                   const double* x,
                                   const double* y,
                                   const double* z,
                                   double* fx,
                                   double* fy,
                                   double* fz) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  double fxi = 0.0;
  double fyi = 0.0;
  double fzi = 0.0;

  double xi = x[i];
  double yi = y[i];
  double zi = z[i];
  double mi = mass[i];

  for (size_t j = 0; j < n; ++j) {
    if (j == i) continue;
    double dx = xi - x[j];
    double dy = yi - y[j];
    double dz = zi - z[j];
    double dist2 = dx*dx + dy*dy + dz*dz + SOFTENING;
    double dist = sqrt(dist2);
    // magnitude of force
    double F = G * mi * mass[j] / dist2;
    // normalize direction safely
    double invnorm = 1.0 / (dist + 1e-30);
    double nx = dx * invnorm;
    double ny = dy * invnorm;
    double nz = dz * invnorm;
    fxi += nx * F;
    fyi += ny * F;
    fzi += nz * F;
  }

  fx[i] = fxi;
  fy[i] = fyi;
  fz[i] = fzi;
}

// integrate kernel to update velocities and positions
__global__ void integrate_gpu(size_t n, double dt,
                              const double* mass,
                              double* x, double* y, double* z,
                              double* vx, double* vy, double* vz,
                              const double* fx, const double* fy, const double* fz) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  double ax = fx[i] / mass[i];
  double ay = fy[i] / mass[i];
  double az = fz[i] / mass[i];

  vx[i] += ax * dt;
  vy[i] += ay * dt;
  vz[i] += az * dt;

  x[i] += vx[i] * dt;
  y[i] += vy[i] * dt;
  z[i] += vz[i] * dt;
}

// helper to run GPU version
double run_gpu(simulation &s, double dt, size_t nbstep, size_t printevery, int blocksize) {
  size_t n = s.nbpart;
  size_t bytes = n * sizeof(double);

  // device pointers
  double *d_mass=nullptr, *d_x=nullptr, *d_y=nullptr, *d_z=nullptr;
  double *d_vx=nullptr, *d_vy=nullptr, *d_vz=nullptr;
  double *d_fx=nullptr, *d_fy=nullptr, *d_fz=nullptr;

  // check n not zero
  if (n == 0) return 0.0;

  CUDA_CHECK(cudaMalloc(&d_mass, bytes));
  CUDA_CHECK(cudaMalloc(&d_x, bytes));
  CUDA_CHECK(cudaMalloc(&d_y, bytes));
  CUDA_CHECK(cudaMalloc(&d_z, bytes));
  CUDA_CHECK(cudaMalloc(&d_vx, bytes));
  CUDA_CHECK(cudaMalloc(&d_vy, bytes));
  CUDA_CHECK(cudaMalloc(&d_vz, bytes));
  CUDA_CHECK(cudaMalloc(&d_fx, bytes));
  CUDA_CHECK(cudaMalloc(&d_fy, bytes));
  CUDA_CHECK(cudaMalloc(&d_fz, bytes));

  // copy initial data to device
  CUDA_CHECK(cudaMemcpy(d_mass, s.mass.data(), bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_x, s.x.data(), bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_y, s.y.data(), bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_z, s.z.data(), bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_vx, s.vx.data(), bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_vy, s.vy.data(), bytes, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_vz, s.vz.data(), bytes, cudaMemcpyHostToDevice));

  // zero forces
  CUDA_CHECK(cudaMemset(d_fx, 0, bytes));
  CUDA_CHECK(cudaMemset(d_fy, 0, bytes));
  CUDA_CHECK(cudaMemset(d_fz, 0, bytes));

  int threads = blocksize > 0 ? blocksize : 128;
  int blocks = (int)((n + threads - 1) / threads);

  // create CUDA events to measure GPU time accurately
  cudaEvent_t startEvent, stopEvent;
  CUDA_CHECK(cudaEventCreate(&startEvent));
  CUDA_CHECK(cudaEventCreate(&stopEvent));
  CUDA_CHECK(cudaEventRecord(startEvent, 0));

  for (size_t step = 0; step < nbstep; ++step) {
    // compute forces
    CUDA_CHECK(cudaMemset(d_fx, 0, bytes));
    CUDA_CHECK(cudaMemset(d_fy, 0, bytes));
    CUDA_CHECK(cudaMemset(d_fz, 0, bytes));

    compute_forces_gpu<<<blocks, threads>>>(n, d_mass, d_x, d_y, d_z, d_fx, d_fy, d_fz);
    // check kernel launch error
    CUDA_CHECK(cudaPeekAtLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // integrate
    integrate_gpu<<<blocks, threads>>>(n, dt, d_mass, d_x, d_y, d_z, d_vx, d_vy, d_vz, d_fx, d_fy, d_fz);
    CUDA_CHECK(cudaPeekAtLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // copy back and print at intervals
    if (step % printevery == 0) {
      CUDA_CHECK(cudaMemcpy(s.x.data(), d_x, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.y.data(), d_y, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.z.data(), d_z, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.vx.data(), d_vx, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.vy.data(), d_vy, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.vz.data(), d_vz, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.fx.data(), d_fx, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.fy.data(), d_fy, bytes, cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(s.fz.data(), d_fz, bytes, cudaMemcpyDeviceToHost));
      dump_state(s);
    }
  }

  CUDA_CHECK(cudaEventRecord(stopEvent, 0));
  CUDA_CHECK(cudaEventSynchronize(stopEvent));
  float ms = 0.0f;
  CUDA_CHECK(cudaEventElapsedTime(&ms, startEvent, stopEvent));

  // cleanup events
  CUDA_CHECK(cudaEventDestroy(startEvent));
  CUDA_CHECK(cudaEventDestroy(stopEvent));

  // free
  CUDA_CHECK(cudaFree(d_mass)); CUDA_CHECK(cudaFree(d_x)); CUDA_CHECK(cudaFree(d_y)); CUDA_CHECK(cudaFree(d_z));
  CUDA_CHECK(cudaFree(d_vx)); CUDA_CHECK(cudaFree(d_vy)); CUDA_CHECK(cudaFree(d_vz));
  CUDA_CHECK(cudaFree(d_fx)); CUDA_CHECK(cudaFree(d_fy)); CUDA_CHECK(cudaFree(d_fz));

  // return seconds measured by events
  return double(ms) * 1e-3;
}

int main(int argc, char* argv[]) {
  if (argc != 5 && argc != 6) {
    std::cerr
      <<"usage: "<<argv[0]<<" <input> <dt> <nbstep> <printevery> [--gpu]"<<"\n"
      <<"input can be:"<<"\n"
      <<"a number (random initialization)"<<"\n"
      <<"planet (initialize with solar system)"<<"\n"
      <<"a filename (load from file in singleline tsv)"<<"\n"
      <<"optional 5th argument: --gpu to run on GPU (requires nvcc build and CUDA-capable device)\n";
    return -1;
  }
  
  double dt = std::atof(argv[2]); //in seconds
  size_t nbstep = std::atol(argv[3]);
  size_t printevery = std::atol(argv[4]);

  // detect optional GPU flag
  bool useGPU = false;
  if (argc == 6) {
    std::string flag = argv[5];
    if (flag == std::string("--gpu")) useGPU = true;
  }
  
  simulation s(1);

  //parse command line
  {
    size_t nbpart = std::atol(argv[1]); //return 0 if not a number
    if ( nbpart > 0) {
      s = simulation(nbpart);
      random_init(s);
    } else {
      std::string inputparam = argv[1];
      if (inputparam == "planet") {
    init_solar(s);
      } else{
    load_from_file(s, inputparam);
      }
    }    
  }

  int blocksize = 128;

  // run either CPU or GPU and time both
  if (useGPU) {
    // print info so logs show we chose GPU path
    std::cerr << "info: running GPU path with blocksize=" << blocksize << " n=" << s.nbpart << "\n";
    double elapsed = run_gpu(s, dt, nbstep, printevery, blocksize);
    std::cerr << "runtime(GPU): " << elapsed << " s\n";
    return 0;
  } else {
    double elapsed = run_cpu(s, dt, nbstep, printevery);
    std::cerr << "runtime(CPU): " << elapsed << " s\n";
    return 0;
  }

  return 0;
}