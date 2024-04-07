#include "stdafx.h"
#include "MathRect.h"
#include "math.h"


CMathRect::CMathRect(double x1, double x2, double y1, double y2, int nx, int ny)
{
	//конструктор с установкой размеров мира
	WXmin = x1; WXmax = x2;
	WYmin = y1; WYmax = y2;
	RM = TM = BM = 20;
	LM = 40;
	RM = 30;

	TM = 15; BM = 15;
	/*
	nx = 7; ny = 7;
	dx = (WXmax - WXmin) / nx;
	dy = (WYmax - WYmin) / ny;
	nx

	dx1 = 3.14159265359;
	dy1 = 1;
	*/
	//dx = 1;
	//dy = 1;

	this->nx = nx;//(WXmax - WXmin) / dx;
	this->ny = ny;// (WYmax - WYmin) / dy;

	dx = x2 / nx;
	dy = (y2 - y1) / ny;
}


CMathRect::~CMathRect()
{
}

void CMathRect::MinMaxX(double x1, double x2)
{
	//Переустановка границ мира по Х
	WXmin = x1; WXmax = x2;
	dx = (WXmax - WXmin) / nx;


}

void CMathRect::MinMaxY(double y1, double y2)
{
	//Переустановка границ мира по Х
	WYmin = y1; WYmax = y2;
	dx = (WYmax - WYmin) / ny;


}
void CMathRect::Clear()
{
	//Очистка области рисования 

	Rectangle(dc, R0.left - 2, R0.top, R0.right, R0.bottom);

	//Rectangle(dc, R.left-1, R.top-1, R.right+1, R.bottom+1);

}




void CMathRect::SetGrades(int inx, int iny)
{
	//Деление мира по Х и Y
	nx = inx;
	ny = iny;
	dx = (WXmax - WXmin) / nx;
	dy = (WYmax - WYmin) / ny;

}
void CMathRect::XGrid(CString xName, CDC* pdc)
{
	CString c;
	double x;
	int k;

	CFont newfont;
	CFont* oldFont;
	newfont.CreateStockObject(ANSI_VAR_FONT);
	oldFont = pdc->SelectObject(&newfont);
	//Вертикальная сетка
	for (int i = 0; i < nx; i++)
	{
		k = 0;
		VertLine(i);
		//MMoveTo(WYmin, WXmin + i*dx);
		c.Format(_T("%.6f"), (double)(WXmin + i*dx));
		SetBkMode(dc, TRANSPARENT);
		x = i*dx;
		if (abs(x) >= 1)
		{
			while (abs(x) >= 1)
			{
				x /= 10;
				k++;

			}

			TextOut(dc, OutX(WXmin + i*dx), OutY(0), c, k + 2);
		}

		else if (abs(x)>1e-7)
		{


			while (abs(x) < 1)
			{
				x *= 10;
				k++;

			}

			TextOut(dc, OutX(WXmin + i*dx), OutY(0), c, k + 3);




		}
		else if (abs(x)<1e-7) TextOut(dc, OutX(WXmin + i*dx), OutY(0), c, 3);



	}
	pdc->SelectObject(oldFont);
	TextOut(dc, OutX(WXmax) + 5, OutY(0) - 8, xName, xName.GetLength());
}

void CMathRect::YGrid(CString yName, CDC* pdc)
{
	CString c;
	double y;
	int k;
	//Вертикальная сетка
	CFont newfont;
	CFont* oldFont;
	newfont.CreateStockObject(ANSI_VAR_FONT);
	oldFont = pdc->SelectObject(&newfont);

	for (int i = 0; i < ny; i++)
	{
		k = 0;
		if (i != 0 && i != ny)HorLine(i);
		//MMoveTo(WXmin, WYmin + i*dy);
		c.Format(_T("%.3f"), (double)(WYmin + dy*i));
		SetBkMode(dc, TRANSPARENT);
		//y = WYmin + dy*i;
		y = _wtof(c);

		while (abs(y) >= 1)
		{
			y /= 10;
			k++;
		}

		while (abs(y) < 1)
		{
			if (y != 0)y *= 10;
			else break;

			k++;
		}

		if ((WYmin + dy*i) >= 0)TextOut(dc, OutX(0) - 38, OutY(WYmin + dy*i) - 8, c, k + 3);
		else TextOut(dc, OutX(0) - 42, OutY(WYmin + dy*i) - 8, c, k + 4);
	}

	pdc->SelectObject(oldFont);
	TextOut(dc, OutX(0) - 30, OutY(WYmax) - 15, yName, yName.GetLength());
}

void CMathRect::HorLine(int k)
{
	//Проведение горизонтальной линии
	MMoveTo(WXmin, WYmin + dy*k);
	MLineTo(WXmax, WYmin + dy*k);
}

void CMathRect::VertLine(int k)
{
	//Ппроведение вертикальной линии 
	MMoveTo(WXmin + k*dx, WYmin);
	MLineTo(WXmin + k*dx, WYmax);
}


void CMathRect::Scale(HDC idc, RECT iR)
{
	//Подключение к контексту и областифы

	dc = idc;
	R0 = iR;
	R.left = R0.left + LM;// +LM;
	R.right = R0.right - RM;// -RM;
	R.top = R0.top + TM;// +TM;
	R.bottom = R0.bottom - BM;// -BM;
	mx = (R.right - R.left) / (WXmax - WXmin);
	my = (R.top - R.bottom) / (WYmax - WYmin);

}
int CMathRect::OutX(double x)
{
	//Пересчет координат по X
	double res = R.left + (int)floor((x - WXmin)*mx);
	return R.left + (int)floor((x - WXmin)*mx);

}

int CMathRect::OutY(double y)
{

	//Пересчет координат по Y
	double res = R.bottom + (int)floor((y - WYmin)*my);
	return R.bottom + (int)floor((y - WYmin)*my);

}

void CMathRect::MMoveTo(double x, double y)
{
	//перенесение пера в мировую точку

	MoveToEx(dc, OutX(x), OutY(y), NULL);

}

void CMathRect::MLineTo(double x, double y)
{
	//проведение линии в мировую точку
	LineTo(dc, OutX(x), OutY(y));


}

void CMathRect::Frame()
{
	//Рамка вокруг рисунка

	//	Rectangle(dc, R.left, R.top, 0, 0);

	//Rectangle(dc, R.left, R.top, R.right, R.bottom);

	//Rectangle(dc, R.left+100, R.top, R.right, R.bottom);

	//Rectangle(dc, R.left, R.top, R.right, R.bottom);
	Rectangle(dc, R.left - 1, R.top - 1, R.right + 1, R.bottom + 1);
}

void CMathRect::XAxis(RECT iR)
{
	//Проведение оси X
	R0 = iR;
	MoveToEx(dc, R.left, OutY(0), NULL);
	LineTo(dc, R.right, OutY(0));
	//LineTo(dc, R.right, (R0.bottom - R0.top) / 2);


}

void CMathRect::YAxis(RECT iR)
{
	//Проведение оси Y(заготовка)
	R0 = iR;

	MoveToEx(dc, OutX(0), R.bottom, NULL);
	LineTo(dc, OutX(0), R.top);

}




void CMathRect::XGrid1()
{
	//Вертикальная сетка(штрих)
	for (int i = 0; i <= nx; i++)
	{

		MMoveTo(WXmin + i*dx + 0.5, WYmin);
		MLineTo(WXmin + i*dx + 0.5, WYmax);
	}
}

void CMathRect::YGrid1()
{
	//Вертикальная сетка(штрих)
	for (int i = 0; i <= ny; i++)
	{
		MMoveTo(WXmin, WYmin + i*dy + 0.5);
		MLineTo(WXmax, WYmin + i*dy + 0.5);
	}
}

double CMathRect::GetBgnPoint()
{
	return WXmin;
}

void CMathRect::XAxis()
{
	MMoveTo(WXmin, 0);
	MLineTo(WXmax, 0);
}

void CMathRect::SetDC(HDC idc)
{
	dc = idc;
}

//void CMathRect::Plot(double dx, double* y, int N)
//{
//	MMoveTo(0, y[0]);
//	for (int i = 1; i < N; i++)
//	{
//		MLineTo(dx*i, y[i]);
//	}
//}


void CMathRect::Plot(double dx, double* y, int N, int shift)
{
	MMoveTo(shift*dx, y[0]);
	if (N == 1 && y[0] >= WYmin && y[0] <= WYmax)
	{
		MLineTo(dx, *y);
	}

	double dx_b;
	int ind;
	for (int i = shift; i < N; i++)
	{
		ind = i - shift;
		if (y[i - 1] <= WYmax&&y[i] > WYmax)
		{

			dx_b = abs((WYmax - y[ind - 1]) / (y[ind] - y[ind - 1])*dx);

			MLineTo(dx*(i - 1) + dx_b, WYmax);
			MMoveTo(dx*i, y[ind]);

		}

		else if (y[ind - 1] > WYmax&&y[ind] > WYmax)
		{

			MMoveTo(dx, y[i]);

		}

		else if (y[ind - 1] > WYmax&&y[ind] <= WYmax)
		{

			dx_b = abs((WYmax - y[i - 1]) / (y[ind] - y[ind - 1])*dx);
			MMoveTo(dx*(i - 1) + dx_b, WYmax);
			MLineTo(dx*i, y[ind]);
		}
		else MLineTo(dx*i, y[ind]);



	}



}

void CMathRect::CreatePlotRect(CDC* dc, CRect& rect_, int left, int top, int right, int bottom, bool bClear, int div, CString xName, CString yName)
{
	//mr = new CMathRect(0, 1 / f_d*(N - 1), -4, 4);

	CRect rect = rect_;
	//GetClientRect(&rect);
	rect.left = rect.left + left;
	rect.top += 140;

	int height = rect.bottom - rect.top;
	rect.bottom = rect.bottom - bottom*((height) / div) - 0;
	rect.top = rect.top + top*((height) / div);
	rect.right = rect.right - right;


	CPen pen(PS_SOLID, 1, ::GetSysColor(COLOR_3DFACE));
	CBrush brush(::GetSysColor(COLOR_3DFACE));

	//CClientDC dc(this);

	this->Scale(*dc, rect);


	CPen* oldPen = dc->SelectObject(&pen);
	CBrush* oldBrush = dc->SelectObject(&brush);

	if (bClear)
	{
		this->Clear();
		dc->SelectObject(oldPen);
		dc->SelectObject(oldBrush);
	}
	pen.DeleteObject();
	brush.DeleteObject();
	brush.CreateSolidBrush(RGB(255, 255, 255));
	pen.CreatePen(PS_SOLID, 1, RGB(0, 0, 0));

	oldPen = dc->SelectObject(&pen);
	oldBrush = dc->SelectObject(&brush);


	if (bClear) this->Frame();

	dc->SelectObject(oldPen);
	dc->SelectObject(oldBrush);

	pen.DeleteObject();
	pen.CreatePen(PS_DOT, 1, RGB(225, 225, 225));
	oldPen = dc->SelectObject(&pen);

	if (bClear)
	{
		this->YGrid(yName, dc);
		this->XGrid(xName, dc);
	}


	dc->SelectObject(oldPen);

	pen.DeleteObject();
	pen.CreatePen(PS_SOLID, 1, RGB(0, 0, 0));
	oldPen = dc->SelectObject(&pen);



	if (bClear) this->XAxis();

}

void CMathRect::DotPlot(double dx, double* y, int N, int shift)
{


	MMoveTo(0, y[0]);
	for (int i = 0; i < N; i++)
	{
		Ellipse(dc, OutX(dx*i)-3, OutY(y[i])-3, OutX(dx*i) + 3, OutY(y[i]) + 3);
	}


}