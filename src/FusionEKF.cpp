#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    /**
    TODO:
    * Finish initializing the FusionEKF.
    */

    float sigma_2_laser = 0.0225;
    float sigma_2_radar_rho = 0.0255;
    float sigma_2_radar_phi = 0.0255;
    float sigma_2_radar_speed = 0.0255;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << sigma_2_laser, 0,
                0, sigma_2_laser;

    R_radar_ = MatrixXd(3, 3);
    R_radar_ << sigma_2_radar_rho, 0, 0,
                0, sigma_2_radar_phi, 0,
                0, 0, sigma_2_radar_speed;


    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    Hj_ = MatrixXd(3, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;
        previous_timestamp_ = measurement_pack.timestamp_;

        //create a 4D state vector, we don't know yet the values of the x_ state
        ekf_.x_ = VectorXd(4);
        //create a 4x4 F_ matrix, still time unspecific
        ekf_.F_ = MatrixXd(4, 4);
        ekf_.F_ <<  1, 0, 1, 0,
                    0, 1, 0, 1,
                    0, 0, 1, 0,
                    0, 0, 0, 1;

        //state covariance matrix P
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ <<  1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1000, 0,
                    0, 0, 0, 1000;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            float rho_dot = measurement_pack.raw_measurements_[2];

            float px = cos(phi) * rho;
            float py = sin(phi) * rho;
            float vx = cos(phi) * rho_dot;
            float vy = sin(phi) * rho_dot;
            ekf_.x_ << px, py, 0, 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            float px = measurement_pack.raw_measurements_[0];
            float py = measurement_pack.raw_measurements_[1];
            float vx = 0;
            float vy = 0;
            ekf_.x_ << px, py, vx, vy;
        }

        // done initializing, no need to predict or update
        if (ekf_.x_[0] != 0 || ekf_.x_[1]!=0){
            is_initialized_ = true;
        }
        return;
    }

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/

    /**
    TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
    */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    //cout << dt << endl;
    previous_timestamp_ = measurement_pack.timestamp_;

    float dt_2 = pow(dt,2);
    float dt_3 = pow(dt,3);
    float dt_4 = pow(dt,4);

    float noise_ax = 5;
    float noise_ay = 5;

    //Modify the F matrix so that the time is integrated
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    //set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  (dt_4/4)*noise_ax, 0, (dt_3/2)*noise_ax, 0,
                0, (dt_4/4)*noise_ay, 0, (dt_3/2)*noise_ay,
                (dt_3/2)*noise_ax, 0, dt_2*noise_ax, 0,
                0, (dt_3/2)*noise_ay, 0, dt_2*noise_ay;

    ekf_.Predict();

    /*****************************************************************************
    *  Update
    ****************************************************************************/

    /**
    TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
    */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        ekf_.H_ = Hj_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
    // Laser updates
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    //cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
}
