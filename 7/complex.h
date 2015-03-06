class complex
{
private:
	double re, im;
public:
	complex();
	complex(double re_, double im_, bool flag);
	void set_re(double re_);
	void set_im(double im_);
	double get_re();
	double get_im();
	void print(bool flag);
	void input();

	complex sum(complex B);
	complex sub(complex B);
	complex mult(complex B);
	complex div(complex B);
	complex conj();
	double sq_mod();
	double modul();
	complex operator + (complex A);
	complex operator * (complex A);
	complex operator / (double);
    void operator += (complex A);
};
