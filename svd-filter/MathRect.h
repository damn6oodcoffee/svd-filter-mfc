#pragma once

class CMathRect
{
	HDC dc; RECT R0, R;
	double WXmin, WXmax, WYmin, WYmax;
	double mx, my;

	int LM, RM, TM, BM;
public:
	int nx, ny;
	double dx, dy;




public:
	BOOL Framed;
	CMathRect(double x1, double x2, double y1, double y2, int nx, int ny);
	virtual ~CMathRect();
	void Scale(HDC idc, RECT iR);
	void MinMaxX(double x1, double x2);
	void MinMaxY(double y1, double y2);
	void Clear();
	int OutX(double x);
	int OutY(double y);
	void MMoveTo(double x, double y);
	void MLineTo(double x, double y);
	void Frame();

	void SetGrades(int inx, int iny);

	void XGrid(CString, CDC*);
	void YGrid(CString, CDC*);


	void XAxis(RECT iR);
	void YAxis(RECT iR);

	void XAxis();

	double GetWX(){ return WXmin; }
	double XWidth()	{ return WXmax - WXmin; }

	void XGrid1();
	void YGrid1();

	double GetBgnPoint();
	void SetDC(HDC idc);

	void CreatePlotRect(CDC* dc, CRect& rect_, int left, int top, int right, int bottom, bool bClear, int div, CString, CString);

	void Plot(double dx, double* y, int N, int shift);
	void DotPlot(double dx, double* y, int N, int shift);
private:
	void HorLine(int k);
	void VertLine(int k);


};

