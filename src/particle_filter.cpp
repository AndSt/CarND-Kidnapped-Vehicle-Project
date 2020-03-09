/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */

   num_particles = 18;  // TODO: Set the number of particles

   std::normal_distribution<double> x_dist(x,std[0]);
   std::normal_distribution<double> y_dist(y,std[1]);
   std::normal_distribution<double> theta_dist(theta,std[2]);

   for (int i = 0; i < num_particles; ++i) {
     Particle p;
     p.id = i;
     p.x = x_dist(gen);
     p.y = y_dist(gen);
     p.theta = theta_dist(gen);
     p.weight = 1;
     particles.push_back(p);
   }
   is_initialized = true;
   std::cout << "Finished init" << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /** transformed_obs
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

   //std::cout << "Start predict" << std::endl;

   std::normal_distribution<double> x_dist(0.0,std_pos[0]);
   std::normal_distribution<double> y_dist(0.0,std_pos[1]);
   std::normal_distribution<double> theta_dist(0.0,std_pos[2]);

   double new_theta = 0;
   for (unsigned int i = 0; i < particles.size(); i++){
     if (fabs(yaw_rate) < 1e-4) {
       particles[i].x += velocity * delta_t * cos(particles[i].theta);
       particles[i].y += velocity * delta_t * sin(particles[i].theta);
     }
     else {
       new_theta = particles[i].theta + yaw_rate * delta_t;
       particles[i].x += (velocity / yaw_rate) * (sin(new_theta) - sin(particles[i].theta));
       particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(new_theta));
       particles[i].theta = new_theta;
     }

     particles[i].x += x_dist(gen);
     particles[i].y += y_dist(gen);
     particles[i].theta += theta_dist(gen);
   }

   //std::cout << "Finished predict" << std::endl;
}

void ParticleFilter::dataAssociation(vector<LandmarkObs>& predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */

   double min_distance = 1e99;
   double distance = 0;

   for(unsigned int i = 0; i < observations.size(); i++){
     for(unsigned int j = 0; j < predicted.size(); j++){
       distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
       if(distance < min_distance){
         min_distance = distance;
         observations[i].id = predicted[j].id;
       }
     }
   }

}

LandmarkObs ParticleFilter::transform_observation(const Particle& p, const LandmarkObs obs){
  double global_x = p.x + obs.x * cos(p.theta) - obs.y * sin(p.theta);
  double global_y = p.y + obs.x * sin(p.theta) + obs.y * cos(p.theta);
  return LandmarkObs{obs.id, global_x, global_y};
}

double ParticleFilter::multiv_prob(double sigma[], LandmarkObs mu, LandmarkObs x) {
  // calculate normalization term
  double gauss_norm = 1.0 / (2 * M_PI * sigma[0] * sigma[1]);

  // calculate exponent
  double exponent = (pow(x.x - mu.x, 2) / (2 * pow(sigma[0], 2))) + (pow(x.y - mu.y, 2) / (2 * pow(sigma[1], 2)));

  // calculate weight using normalization terms and exponent
  double weight = gauss_norm * exp(-exponent);
  return weight;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

   // Algorithm:
   // For each particle:
   //   Find suting LandmarkObs
   //   Align coordinate systems
   //   compute Association
   //   update using Multivariate_normal_distribution

   for (unsigned int i = 0; i < particles.size(); i++){

     Particle p = particles[i];
     vector<LandmarkObs> predicted_landmarks;

     for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++){
       LandmarkObs lmo = LandmarkObs();
       lmo.id = map_landmarks.landmark_list[j].id_i;
       lmo.x = map_landmarks.landmark_list[j].x_f;
       lmo.y = map_landmarks.landmark_list[j].y_f;
       double lm_dist = dist(lmo.x, lmo.y, p.x, p.y);
       if (lm_dist <= sensor_range) {
         predicted_landmarks.push_back(lmo);
       }
     }

     vector<LandmarkObs> transformed_obs;
     for (int j = 0; j < observations.size(); j++){
       transformed_obs.push_back(transform_observation(p, observations[j]));
     }

     dataAssociation(predicted_landmarks, transformed_obs);

     double weight = 1.0;

     //std::cout << "num_trafos:" << transformed_obs.size() << std::endl;
     //std::cout << "num_lms:" << predicted_landmarks.size() << std::endl;
     for (int j=0; j< transformed_obs.size(); j++){
       for (int k=0; k < predicted_landmarks.size(); k++){
         if (predicted_landmarks[k].id == transformed_obs[j].id) {
           weight *= multiv_prob(std_landmark, predicted_landmarks[k], transformed_obs[j]);
           break;
         }
       }
     }
     particles[i].weight = weight;
   }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   vector<Particle> resampled_particles;
   vector<double> weights;
   for(int i=0; i<particles.size();i++){
       weights.push_back(particles[i].weight);
   }

   std::random_device rd;
   std::default_random_engine gen(rd());

   std::uniform_int_distribution<int> ind_dist(0, particles.size()-1);
   int index = ind_dist(gen);

   double max_weight = *max_element(weights.begin(), weights.end());
   std::uniform_real_distribution<double> beta_dist(0.0, max_weight);

   double beta = 0.0;

   for (int i = 0; i < particles.size(); i++) {
       beta += beta_dist(gen) * 2.0;
       while (beta > weights[index]) {
           beta -= weights[index];
           index = (index + 1) % particles.size();
       }
       resampled_particles.push_back(particles[index]);
   }

   particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
