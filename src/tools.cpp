#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    // ... your code here
    if (estimations.size() == 0)
    {
        cout << "estimation vector size is zero";
    }else if(estimations.size() != ground_truth.size())
    {
        cout << "estimation vector size is not equal to ground truth vector size";
    }

    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd est = estimations[i];
        VectorXd gt = ground_truth[i];

        VectorXd a = est - gt;
        VectorXd b = a.array() * a.array();

        rmse += b;
    }

    //calculate the mean
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    float val1 = pow(px,2) + pow(py,2);
    float val2 = sqrt(val1);
    float val3 = (val1*val2);

    //check division by zero
    if (val1 == 0)
    {
        cout << "CalculateJacobian () - Error - Division by Zero";
    }else{
        //compute the Jacobian matrix
        Hj(0,0) = px / val2;
        Hj(0,1) = py / val2;
        Hj(1,0) = -py / val1;
        Hj(1,1) = px / val1;
        Hj(2,0) = py*(vx*py - vy*px) / val3;
        Hj(2,1) = px*(vy*px - vx*py) / val3;
        Hj(2,2) = px / val2;
        Hj(2,3) = py / val2;
    }

    return Hj;
}
