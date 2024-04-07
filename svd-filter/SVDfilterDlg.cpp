
// SVDfilterDlg.cpp : implementation file
//

#include "stdafx.h"
#include "SVDfilter.h"
#include "SVDfilterDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CSVDfilterDlg dialog



CSVDfilterDlg::CSVDfilterDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CSVDfilterDlg::IDD, pParent)
	, m_N(1024)
	, m_fd(120)
	, m_L(100)
	, m_a1(1.7)
	, m_a2(1.3)
	, m_a3(1.4)
	, m_f1(11.2)
	, m_f2(33.7)
	, m_f3(48.7)
	, m_phi1(22)
	, m_phi2(12)
	, m_phi3(1)
	, m_noise(0)
	, m_thrhold(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CSVDfilterDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_N, m_N);
	DDX_Text(pDX, IDC_FD, m_fd);
	DDX_Text(pDX, IDC_L, m_L);
	DDX_Text(pDX, IDC_A1, m_a1);
	DDX_Text(pDX, IDC_A2, m_a2);
	DDX_Text(pDX, IDC_A3, m_a3);
	DDX_Text(pDX, IDC_F1, m_f1);
	DDX_Text(pDX, IDC_F2, m_f2);
	DDX_Text(pDX, IDC_F3, m_f3);
	DDX_Text(pDX, IDC_PHI1, m_phi1);
	DDX_Text(pDX, IDC_PHI2, m_phi2);
	DDX_Text(pDX, IDC_PHI3, m_phi3);
	DDX_Text(pDX, IDC_NOISECOEF, m_noise);
	DDX_Text(pDX, IDC_THRHOLD, m_thrhold);
	DDX_Control(pDX, IDC_NOISECOEF, m_noiseEditBox);
	DDX_Control(pDX, IDC_THRHOLD, m_thrEditBox);
}

BEGIN_MESSAGE_MAP(CSVDfilterDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_RUN, &CSVDfilterDlg::OnBnClickedRun)
	ON_BN_CLICKED(IDC_FILTER, &CSVDfilterDlg::OnBnClickedFilter)
	ON_BN_CLICKED(IDC_CHECK1, &CSVDfilterDlg::OnBnClickedCheck1)
END_MESSAGE_MAP()


// CSVDfilterDlg message handlers

BOOL CSVDfilterDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

	dc = new CClientDC(this);
	memDC = new CDC();

	memDCbuf = new CDC();

	CRect rect;
	GetClientRect(&rect);
	
	InitHeight = rect.Height();
	InitWidth = rect.Width();
	m_noiseEditBox.SetReadOnly();
	//m_thrEditBox.SetReadOnly();
	



	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CSVDfilterDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();

		CRect rect;
		GetClientRect(&rect);


		CPen pen(PS_SOLID, 1, RGB(255, 9, 90));
		CPen* oldPen;


		if (bSignalCreated)
		{

		//CPen pen(PS_SOLID, 1, RGB(215, 117, 234));
		//CPen* oldPen;
		//oldPen = memDC->SelectObject(&pen);
		//memDC->BitBlt(0, 0, rect.Width(), rect.Height(), memDCbuf, 0, 0, SRCCOPY);
		//
		//mr1->MMoveTo(0, s[0]);
		//mr1->Plot(1 / m_fd, s, m_N, 0);
			memDC->BitBlt(0, 0, rect.Width(), rect.Height(), memDCbuf, 0, 0, SRCCOPY);

		
			if (bFiltered){
				CPen *oldPen;
				CPen pen2(PS_SOLID, 2, RGB(255, 71, 255));

				mr2->MMoveTo(0, P[0]);
				oldPen = memDC->SelectObject(&pen2);
				mr2->Plot(m_fd / m_N, P, m_N / 2 + 1, 0);
				memDC->SelectObject(oldPen);

			}
			CBrush brush1(RGB(0, 0, 250));
			CPen pen3(PS_SOLID, 2, RGB(255, 0, 0));
		
	
		//	mr3 = new CMathRect(0, m_L, 0, max_a + max_a / 20, 16, 8);
		//	mr3->CreatePlotRect(memDC, rect, 0, 2, 0, 0, true, 3, _T("n"), _T("sigma(n)"));
			mr3->MMoveTo(0, sigma[0]);

		

			CBrush* pOldbrush = memDC->SelectObject(&brush1);
			oldPen = memDC->SelectObject(&pen3);
			mr3->DotPlot(1, sigma, m_L, 0);

			memDC->SelectObject(oldPen);
			memDC->SelectObject(pOldbrush);


			dc->StretchBlt(0, 140, rect.Width(), rect.Height() - 140, memDC, 0, 140, InitWidth, InitHeight - 140, SRCCOPY);
		//	memDC->SelectObject(oldPen);

		}


	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CSVDfilterDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CSVDfilterDlg::OnBnClickedRun()
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE);

	if (bSignalCreated)
	{
		delete[] s;
		delete[] r_s;
		delete[] P;

		delete[] u;
		delete[] v;
		delete[] sigma;
		delete[] r_s_filtered;
	}
	s = new double[m_N];

	r_s_filtered = new double[m_L*m_L];

	r_s = new double[m_L * m_L];
	P = new double[m_N];

	u = new double[m_L*m_L];
	v = new double[m_L*m_L];
	sigma = new double[m_L];
	sigma_buf = new double[m_L];

	SignalParam p;
	p.a1 = m_a1; p.f1 = m_f1; p.phi1 = m_phi1;
	
	p.a2 = m_a2; p.f2 = m_f2; p.phi2 = m_phi2;

	p.a3 = m_a3; p.f3 = m_f3; p.phi3 = m_phi3;

	bSignalCreated = GenerateSineSignal(s, m_N, m_fd, p);
	

	if(bAddNoise)AddNoise(s, m_N, m_noise);

	Autocorr(s, r_s, m_N, m_L);
	CorrelogramPSD(r_s, P, m_N , m_L, m_fd);

	svd_double(m_L, m_L, r_s, u, v, sigma);

	for (int i = 0; i < m_L; i++){ sigma_buf[i] = sigma[i]; }
	//PrintSigma(sigma, m_L, m_L, u, v);
	//Filter(sigma, u, v, r_s_filtered, m_thrhold, m_L);


	CRect rect;
	GetClientRect(&rect);

	delete memDC;
	memDC = new CDC();

	delete memDCbuf;
	memDCbuf = new CDC();

	CBitmap bitmap;
	bitmap.CreateCompatibleBitmap(dc, rect.Width(), rect.Height());

	CBitmap bitmap1;
	bitmap1.CreateCompatibleBitmap(dc, rect.Width(), rect.Height());
	
	memDC->CreateCompatibleDC(dc);
	memDCbuf->CreateCompatibleDC(dc);

	CBrush brush(::GetSysColor(COLOR_3DFACE));



	CPen pen1(PS_SOLID, 1, RGB(255, 0, 90));
	//CPen pen2(PS_SOLID, 2, RGB(0, 0, 255));
	CPen pen2(PS_SOLID, 2, RGB(0, 0, 0));

	CBitmap* pOldBitmap = memDC->SelectObject(&bitmap);
	CBitmap* pOldBitmap1 = memDCbuf->SelectObject(&bitmap1);


	CPen *oldPen;





	double max_a = abs(s[0]);
	for (int i = 0; i < m_N; i++)
	{
		if (max_a < abs(s[i])) max_a = abs(s[i]);
	}

	
	mr1 = new CMathRect(0, (m_N - 1) / m_fd, -(max_a + max_a / 20), max_a + max_a / 20, 16, 8);
	mr1->CreatePlotRect(memDC, rect, 0, 0, 0, 2, true, 3, _T("t"), _T("s(t)"));
	mr1->MMoveTo(0, s[0]);

	oldPen = memDC->SelectObject(&pen1);
	mr1->Plot(1 / m_fd, s, m_N, 0);

	memDC->SelectObject(oldPen);
	

	 max_a = P[0];
	for (int i = 0; i < m_N; i++)
	{
		if (max_a < P[i]) max_a = P[i];
	}

	mr2 = new CMathRect(0, m_fd/2, -(max_a + max_a / 20), max_a + max_a / 20, 16, 8);
	mr2->CreatePlotRect(memDC, rect, 0,1, 0,1, true,3, _T("f"), _T("P(f)"));
	mr2->MMoveTo(0, P[0]);

	oldPen = memDC->SelectObject(&pen2);
	mr2->Plot(m_fd/m_N, P, m_N/2+1, 0);

	memDC->SelectObject(oldPen);


	mr3 = new CMathRect(0, m_L, 0, sigma[0] + sigma[0] / 20, 16, 8);
	mr3->CreatePlotRect(memDC, rect, 0, 2, 0, 0, true, 3, _T("n"), _T("sigma(n)"));


	CPen penLine1(PS_SOLID, 2, RGB(0, 0, 0));
	CPen penLine2(PS_SOLID, 2, RGB(255, 71, 255));
	oldPen = memDC->SelectObject(&penLine1);
	memDC->MoveTo(800, 350);
	memDC->LineTo(850, 350);
	memDC->TextOutW(853, 340, L"ÑÏÌ ñ øóìîì");
	memDC->SelectObject(&penLine2);
	memDC->MoveTo(950, 350);
	memDC->LineTo(1000, 350);
	memDC->TextOutW(1003, 340, L"Îòôèëüòð. ÑÏÌ");
	memDC->SelectObject(oldPen);

	memDCbuf->BitBlt(0, 0, rect.Width(), rect.Height(), memDC, 0, 0, SRCCOPY);

	
	
	

	OnPaint();
	

}


void CSVDfilterDlg::OnBnClickedFilter()
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE);

	for (int i = 0; i < m_L; i++){ sigma[i] = sigma_buf[i]; }
	Filter(sigma, u, v, r_s_filtered, m_thrhold, m_L);

	CorrelogramPSD(r_s_filtered, P, m_N, m_L, m_fd);


	bFiltered = true;
	OnPaint();
	bFiltered = false;







	dc->MoveTo(mr3->OutX(0), mr3->OutY(sigma[0] * m_thrhold / 100.));
	dc->LineTo(mr3->OutX(m_L), mr3->OutY(sigma[0] * m_thrhold / 100.));




}


void CSVDfilterDlg::OnBnClickedCheck1()
{
	// TODO: Add your control notification handler code here
	bAddNoise = !bAddNoise;
	m_noiseEditBox.SetReadOnly(!bAddNoise);
//	m_thrEditBox.SetReadOnly(!bAddNoise);
}
