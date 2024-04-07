
// SVDfilterDlg.h : header file
//

#pragma once
#include "MathRect.h"
#include "signal.h"
#include "time.h"
#include "afxwin.h"

// CSVDfilterDlg dialog
class CSVDfilterDlg : public CDialogEx
{
// Construction
public:
	CSVDfilterDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_SVDFILTER_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedRun();

	CMathRect* mr1, *mr2, *mr3;
	CClientDC* dc;
	CDC* memDC, *memDCbuf;

	double* s, *r_s, *P;

	double* u, *v, *sigma;
	double* sigma_buf;

	double* r_s_filtered;

	bool bSignalCreated = false;
	bool bFiltered = false;
	bool bAddNoise = false;

	int InitWidth, InitHeight;
	int m_N;
	double m_fd;
	int m_L;
	double m_a1;
	double m_a2;
	double m_a3;
	double m_f1;
	double m_f2;
	double m_f3;
	double m_phi1;
	double m_phi2;
	double m_phi3;
	double m_noise;
	double m_thrhold;
	afx_msg void OnBnClickedFilter();
	afx_msg void OnBnClickedCheck1();
	CEdit m_noiseEditBox;
	CEdit m_thrEditBox;
};
