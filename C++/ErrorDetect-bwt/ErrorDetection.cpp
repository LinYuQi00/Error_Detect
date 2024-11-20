//
// Created by liufei on 2015/12/8.
//

#include <fstream>
#include <iostream>
#include <cmath>
#include <armadillo>
#include "ErrorDetection.h"
#include <string.h>
#include <algorithm>
#include "results.h"

using namespace std;
using namespace arma;

#define METHOD_ID 5 //1:Lpvs  2:LD3  3:LD2  4:LOH  5:TanKeqin
#define OUTPUT_FORM 2 //1:Origin  2:MyResults65  3:MyResults4m0

int ErrorDetection::process(std::string ioe_fileName, std::string gcp_fileName,
                            std::string pixel_fileName, std::string output_fileName){

    /* _____ 1.read data _____ */
    ReadTestFile(ioe_fileName, gcp_fileName, pixel_fileName);

    /* _____ 2.get common pixel point _____ */
    GetCommonPixel();

    /* _____ 3.relative orientation and error detection _____ */
    DetectError(output_fileName);

    return 0;
}

int ErrorDetection::ReadTestFile(std::string ioe_fileName,
                                 std::string gcp_fileName, std::string pixel_fileName) {
    /*_____open ioe_file and check if opened_____*/
    ifstream in_file(ioe_fileName);
    if(!in_file.is_open()){
        cout << "Open inner orientation elements file failed!\n";
        return -1;
    }
    /*_____get IOE_____*/
    string str;
    getline(in_file, str); // read the first line: 1000  4
    Inner in;
    in_file >> in.f >> in.x0 >> in.y0;
    in_file.close();


    /*_____open gcp_file and check if opened_____*/
    in_file.open(gcp_fileName);
    if(!in_file.is_open()){
        cout << "Open ground control point file failed!\n";
        return -1;
    }
    /*_____get data_____*/
    while(true){
        int id{0}, d{0};
        double x{0}, y{0}, z{0};
        in_file >> id >> y >> x >> z ;
        // in_file >> id >> y >> x >> z >> d;
        if(id == -99)
            break;
        GCP point(id, x, y, z);
        //GCP point(id, x, y, z, d);
        GCPs.push_back(point);
    }
    in_file.close();


    /*_____open pixel_file_____*/
    in_file.open(pixel_fileName);
    if(!in_file.is_open()){
        cout << "Open pixel file failed!\n";
        return -1;
    }

    /*_____read pixel data_____*/

    //read first img
    getline(in_file, str);
    int id1{0};
    in_file >> id1;
    if(id1 < 0) id1 = -id1;
    getline(in_file, str);
    while(true){
        int ID{0};
        double x{0}, y{0};
        in_file >> ID >> x >> y;
        ImagePixel pixel{ID, x, y};
        if(ID == -99)
            break;
        img1.insert(pixel);
    }

    //read second img
    getline(in_file, str);
    int id2{0};
    in_file >> id2;
    if(id2 < 0) id2 = -id2;
    getline(in_file, str);
    while(true){
        int ID{0};
        double x{0}, y{0};
        in_file >> ID >> x >> y;
        ImagePixel pixel{ID, x, y};
        if(ID == -99)
            break;
        img2.insert(pixel);
    }
    in_file.close();

    /* _____ Update img information _____ */
    img1.setIOE(in);
    img2.setIOE(in);
    img1.setID(id1);
    img2.setID(id2);

    return 0;
}


int ErrorDetection::GetCommonPixel(){
    /* _____ search img1 pixel in img2 _____ */
    for(auto it = img1.begin(); it != img1.end(); ++it){
        if(img2.find(it->first))
            cPixels.push_back(it->first);
    }

    return 0;
}

int ErrorDetection::DetectError(std::string output_fileName){
    /* _____ init value _____ */
    int num = static_cast<int>(cPixels.size()); // 
    double f = img1.getInner().f;
    double phi{0}, omega{0}, kappa{0}, u{0}, v{0};
    double a1{0}, a2{0}, a3{0}, b1{0}, b2{0}, b3{0}, c1{0}, c2{0}, c3{0};
    mat A(num, 5, fill::zeros);
    mat DX(5, 1, fill::zeros);
    mat P(num, num, fill::zeros);
    mat Q(num, 1, fill::zeros);
    for(int i = 0; i < num; ++i)
        P(i, i) = 1;

    int loopNum1{0};
    /* _____ begin select weight iteration loop _____ */
    while(true){
        ++loopNum1;
        int loopNum2{0};
        /* _____ begin relative orientation _____ */
        while(true){
	    ++loopNum2;

	    /* _____ Bx, By, Bz _____ */
            double Bx = 1;  //set Bx equal 1
            double By = Bx*u;
            double Bz = Bx*v;
            a1 = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa);
            a2 =-cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa);
            a3 =-sin(phi)*cos(omega);
            b1 = cos(omega)*sin(kappa);
            b2 = cos(omega)*cos(kappa);
            b3 =-sin(omega);
            c1 = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa);
            c2 =-sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa);
            c3 = cos(phi)*cos(omega);

            /* _____ one equation for one point _____ */
            for(int i = 0; i < num; ++i){
                /* _____ transform from x,y to X,Y,Z _____ */
                double X1 = img1[cPixels[i]].x;
                double Y1 = img1[cPixels[i]].y;
                double Z1 = -f;
                double x2 = img2[cPixels[i]].x;
                double y2 = img2[cPixels[i]].y;
                double X2 = a1*x2 + a2*y2 - a3*f;
                double Y2 = b1*x2 + b2*y2 - b3*f;
                double Z2 = c1*x2 + c2*y2 - c3*f;

                /* _____ Point projection coefficient _____ */
                double N1 = (Bx*Z2-Bz*X2) / (X1*Z2-Z1*X2);
                double N2 = (Bx*Z1-Bz*X1) / (X1*Z2-Z1*X2);

                /* _____ Matrix A _____ */
                A(i, 0) = -X2*Y2*N2/Z2;
                A(i, 1) = -(Z2+Y2*Y2/Z2)*N2;
                A(i, 2) = X2*N2;
                A(i, 3) = Bx;
                A(i, 4) = -Y2*Bx/Z2;

                /* _____ Matrix Q _____ */
                Q(i, 0) = N1*Y1 - N2*Y2 - By;
            }

            /* _____ Solve equation and corrections _____ */
            DX = (A.t()*P*A).i()*A.t()*P*Q;
            phi += DX(0, 0);
            omega += DX(1, 0);
            kappa += DX(2, 0);
            u += DX(3, 0);
            v += DX(4, 0);

            mat T = abs(DX);
            if(T.max() < 3e-5 || loopNum2 > 30)
                break;
        }

        /* _____ caculate V,Qvv,sigma0 _____ */
        mat V = A*DX - Q;
        mat Qvv = P.i() - A*(A.t()*P*A).i()*A.t();
        mat VtPV = V.t()*P*V;
		double sigma0 = sqrt(VtPV(0, 0) / (num - 5));
		double d = 3.5 + 82 / (81 + pow((sigma0 / 0.0028), 4));

        /* _____ recaculate matrix P _____ */
		mat P_pre = P; float K;
        for(int i = 0; i < num; ++i){
			// F ����ļ���ͳ����T,Ȩ����P
            double Ti = V(i, 0)*V(i, 0) / (sigma0*sigma0*Qvv(i, i)*P(i, i));

   //         if(Ti >=1 && loopNum1 <= 3)
   //             P(i, i) = 1/Ti;
   //         else if(Ti >=3.29*3.29 && loopNum1 >3){
   //             P(i, i) = 1/Ti;
   //             if(P(i, i) < 5e-5)
   //                 ePoints.insert(cPixels[i]);
   //         }


			double wi = fabs(V(i, 0)) / (sigma0 * sqrt(Qvv(i, i)));

			if (loopNum1 <= 3)
				K = 1;
			else
				K = 3.29;
#if METHOD_ID==1
			//ѡȨ��������Ȩ����Lpvs
			if (Ti <= K)
				P(i, i) = 1;
			else
				P(i, i) = 1 / Ti;
			
			if (P(i, i) < 0.001)
				ePoints.insert(cPixels[i]);
#endif

#if METHOD_ID==2
			//�Ľ���ĵ���LD3
			if (Ti <= K)
				P(i, i) = 1;
			else if ((loopNum1 == 2 || loopNum1 == 3) && Ti > K)
				P(i, i) = pow(exp(-pow(Ti, 2.2)), 0.05);
			else if(loopNum1 > 3 && Ti > K)
				P(i, i) = pow(exp(-pow(Ti, 1.5)), 0.05);

			if (P(i, i) < 0.01)
				ePoints.insert(cPixels[i]);
#endif

#if METHOD_ID==3
			//����LD2
			if (loopNum1 == 1 && wi <= 2)
				P(i, i) = 1;
			else if ((loopNum1 == 2 || loopNum1 == 3) && wi > 2)
				P(i, i) = pow(exp(-pow(wi, 4.4)), 0.05);
			else if (loopNum1 > 3 && wi > 2)
				P(i, i) = pow(exp(-pow(wi, 3.0)), 0.05);

			if (P(i, i) < 0.01)
				ePoints.insert(cPixels[i]);
#endif
#if METHOD_ID==4
			//�ݺ���LOH
			if (fabs(wi) <= 2)
				P(i, i) = 1;
			else
				P(i, i) = 1 / wi*wi;

			if (P(i, i) < 0.01)
				ePoints.insert(cPixels[i]);
#endif
#if METHOD_ID==5
			//��Robustԭ��������ѡȨ������
			double Pi_pre;
			if (i > 0)
				Pi_pre = P(i - 1, i - 1);
			else
				Pi_pre = 1;
			double  alpha_i = sqrt(Pi_pre) / (1.4*sigma0*sqrt(Qvv(i, i)*P(i, i)));
			P(i, i) = Pi_pre / (1 + pow((alpha_i*fabs(V(i, 0))),d));


			if (P(i, i) < 0.01)
				ePoints.insert(cPixels[i]);
#endif
        }

        /* _____ break condition _____ */
        mat DP = abs(P_pre - P);
        double max = DP.max();
//#if METHOD_ID==2
//		if (loopNum1 == 1)
//			continue;
//#endif
        if(DP.max() < 0.001 || loopNum2 > 30)
            break;
    }

	
	//double *errorVal = new double[ePoints.size()];
	vector<results> vect;
	results res("",0,0);

	int eNum = 0;

    /* _____ calculate error in y through [phi,omega,kappa,u,v] _____ */
    double Bx = 1;
    double By = Bx*u;
    double Bz = Bx*v;
//    int eNum = static_cast<int>(ePoints.size());
    for(auto it = ePoints.begin(); it != ePoints.end(); ++it){
        double X1 = img1[*it].x;
        double Y1 = img1[*it].y;
        double Z1 = -f;
        double x2 = img2[*it].x;
        double y2 = img2[*it].y;
        double X2 = a1*x2 + a2*y2 - a3*f;
        double Y2 = b1*x2 + b2*y2 - b3*f;
        double Z2 = c1*x2 + c2*y2 - c3*f;

        double N1 = (Bx*Z2-Bz*X2) / (X1*Z2-Z1*X2);
        double N2 = (Bx*Z1-Bz*X1) / (X1*Z2-Z1*X2);

        /* _____ get dy _____ */
        double dy = (N1*Y1-By)/N2-Y2;

		stringstream ss;
		string pid;
		ss << *it;
		ss >> pid;
		res.PointID = pid;
		res.errAbs = fabs(dy / m0);
		res.errVal = dy / m0;
		vect.push_back(res);

		eNum++;
    }
	/* _____ write result in the file  _____ */
	ofstream out_file(output_fileName);
	if (!out_file.is_open()){
		cout << "Create result file failed!\n";
		return -1;
	}

#if OUTPUT_FORM==1
	out_file << "Error Detection result:\n" << "Total error points: " << ePoints.size() << "\n"
		<< "Point ID\t" << "absolute value\t" << "relative to m0\n";
	for (int i = 0; i < vect.size(); i++){
		out_file.setf(ios::left);
		out_file.width(9);
		out_file << vect[i].PointID << "\t";
		out_file.setf(ios::left);
		out_file.width(9);
		out_file << vect[i].errVal*m0 << "\t";
		out_file.setf(ios::left);
		out_file.width(9);
		out_file << vect[i].errVal << "\n";
	}
#endif
#if OUTPUT_FORM==2
	/* sort the value of ePoints */
	stable_sort(vect.begin(), vect.end(), greater<results>());
	/* get the first 65 error points */
	vector<results> newVect;
	for (int i = 0; i < vect.size(); i++){
		if (i < 65){
			res.PointID = vect[i].PointID;
			res.errVal = vect[i].errVal;
			res.errAbs = vect[i].errAbs;
			newVect.push_back(res);
		}
		else{
			break;
		}
	}
	/* sort results by point ID */
	stable_sort(newVect.begin(), newVect.end(), less<results>());
	
	out_file << "Error Detection result:\n" << "65 of Most biggest error points: " << "\n"
		<< "Point ID\t" << "误差大小\t" << "误差倍数\n";
	for (int i = 0; i < newVect.size(); i++){
		out_file.setf(ios::left);
		out_file.width(9);
		out_file << newVect[i].PointID << "\t" << newVect[i].errVal*m0 << "\t" << newVect[i].errVal << "\n";
	}

#endif

#if OUTPUT_FORM==3
	out_file << "Error Detection result:\n" << "Total error points: " << ePoints.size() << "\n"
		<< "Point ID\t" << "absolute value\t" << "relative to m0\n";
	for (int i = 0; i < vect.size(); i++){
		if (vect[i].errAbs >= 4){
			out_file.setf(ios::left);
			out_file.width(9);
			out_file << vect[i].PointID << "\t";
			out_file.setf(ios::left);
			out_file.width(9);
			out_file << vect[i].errVal*m0 << "\t";
			out_file.setf(ios::left);
			out_file.width(9);
			out_file << vect[i].errVal << "\n";
		}
	}
#endif

	out_file.close();
    return 0;
}
