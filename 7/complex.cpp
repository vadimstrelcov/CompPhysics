#include <iostream>
#include "complex.h"
#include <cmath>
#include <cstdlib>
using namespace std;

complex::complex() {
    this->re=0;
    this->im=0;
}

complex::complex(double re_, double im_, bool flag) {
    if (flag) {
        this->re=re_;
        this->im=im_;
    } else {
        this->re=re_*cos(im_);
        this->im=re_*sin(im_);
    }
}

void complex::set_re(double re_) {
	this->re=re_;
}

void complex::set_im(double im_) {
	this->im=im_;
}

double complex::get_re() {
	return this->re;
}

double complex::get_im() {
	return this->im;
}

void complex::print(bool flag) {
    if (flag) {
        cout<<"("<<this->re<<";"<<this->im<<")"<<endl;
    } else {
        cout<<this->modul()<<"exp("<<atan2(this->im,this->re)<<"i)"<<endl;
    }
}

void complex::input() {
	cin>>this->re>>this->im;
}

complex complex::sum(complex B) {
	return complex(this->re+B.get_re(),this->im+B.get_im(),1);
}

complex complex::sub(complex B) {
	return complex(this->re-B.get_re(),this->im-B.get_im(),1);
}

complex complex::mult(complex B) {
	return complex(this->re*B.get_re()-this->im*B.get_im(),this->re*B.get_im()+this->im*B.get_re(),1);
}

complex complex::div(complex B) {
	if (fabs(B.get_re())<1e-10 && fabs(B.get_im())<1e-10) {
		cerr<<"division by 0";
		exit(1);
	} else {
		return complex( (this->re*B.get_re()+this->im*B.get_im())/(B.get_re()*B.get_re()+B.get_im()*B.get_im()),
                 (B.get_re()*this->im-B.get_im()*this->re)/(B.get_re()*B.get_re()+B.get_im()*B.get_im()) ,1);
	}

}

complex complex::conj() {
	return complex(this->re,-this->im,1);
}

double complex::sq_mod() {
    return (this->re*this->re+this->im*this->im);
}

double complex::modul() {
	return sqrt(this->re*this->re+this->im*this->im);
}

complex complex::operator + (complex A) {
	return complex(this->re+A.get_re(),this->im+A.get_im(),1);
}

complex complex::operator * (complex A) {
   	return complex(this->re*A.get_re()-this->im*A.get_im(),this->re*A.get_im()+this->im*A.get_re(),1);
}

complex complex::operator / (double x) {
	if (fabs(x)<1e-10) {
		cerr<<"division by 0";
		exit(1);
	} else {
		return complex(this->re/x,this->im/x,1);
	}
}

void complex::operator += (complex A) {
    this->re+=A.get_re();
    this->im+=A.get_im();
}
