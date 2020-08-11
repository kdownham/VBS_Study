//This is a ROOT macro designed to plot the invariant mass distribution for electroweak WW vector boson scattering and the prominent backgrounds

#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TRint.h"

TCanvas* example_plot( int iPeriod, int iPos );
 
void PlotMacro()
{
   //	gROOT->LoadMacro("tdrstyle.C");
   setTDRStyle();

   //	gROOT->LoadMacro("CMS_lumi.C");
   
   writeExtraText = true;
   extraText = "Preliminary";
   lumi_8TeV = "19.1 fb^{-1}";
   lumi_7TeV = "4.9 fb^{-1}";
   lumi_13TeV = "35.9 fb^{-1}";
   lumi_sqrtS = "13 TeV";

   int iPeriod = 3;

   example_plot( iPeriod, 0 );

}

TCanvas* example_plot( int iPeriod, int iPos){
	
	int W = 800;
	int H = 600;

	int H_ref = 600;
	int W_ref = 800;
	
	//references for T, B, L, R
	float T = 0.08*H_ref;
	float B = 0.12*H_ref; 
        float L = 0.12*W_ref;
        float R = 0.04*W_ref;

        TString canvName = "FigExample_";
  	canvName += W;
  	canvName += "-";
  	canvName += H;
  	canvName += "_";  
  	canvName += iPeriod;
  	if( writeExtraText ) canvName += "-prelim";
  	if( iPos%10==0 ) canvName += "-out";
  	else if( iPos%10==1 ) canvName += "-left";
  	else if( iPos%10==2 )  canvName += "-center";
  	else if( iPos%10==3 )  canvName += "-right";
	
  	TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  	canv->SetFillColor(0);
  	canv->SetBorderMode(0);
  	canv->SetFrameFillStyle(0);
  	canv->SetFrameBorderMode(0);
  	canv->SetLeftMargin( L/W );
  	canv->SetRightMargin( R/W );
  	canv->SetTopMargin( T/H );
  	canv->SetBottomMargin( B/H );
  	canv->SetTickx(0);
  	canv->SetTicky(0);
	
  	TH1* h = new TH1D("h","h",6,500,3500);
  	h->GetXaxis()->SetNdivisions(6,5,0);
  	h->GetXaxis()->SetTitle("m_{jj} (GeV)");  
  	h->GetYaxis()->SetNdivisions(6,5,0);
  	h->GetYaxis()->SetTitleOffset(1);
  	h->GetYaxis()->SetTitle("d#sigma/d#it{m_{jj}} [fb / 500 GeV]");  
	
  	h->SetMaximum( 50 );
  	if( iPos==1 ) h->SetMaximum( 50 );
  	h->Draw();

  	int histLineColor = kOrange+7;
  	int histFillColor = kOrange-2;
  	float markerSize  = 1.0;

	{
    TLatex latex;
				
    int n_ = 5;

    float x1_l = 0.92;
    float y1_l = 0.60;

    float dx_l = 0.30;
    float dy_l = 0.18;
    float x0_l = x1_l-dx_l;
    float y0_l = y1_l-dy_l;

    TPad* legend = new TPad("legend_0","legend_0",x0_l,y0_l,x1_l, y1_l );

    legend->Draw();
    legend->cd();
		
    float ar_l = dy_l/dx_l;
		
    float x_l[1];
    float ex_l[1];
    float y_l[1];
    float ey_l[1];

    float gap_ = 1./(n_+1);
		
    float bwx_ = 0.12;
    float bwy_ = gap_/1.5;
		
    x_l[0] = 1.2*bwx_;

    y_l[0] = 1-gap_;
    ex_l[0] = 0;
    ey_l[0] = 0.04/ar_l;
		
    TGraph* gr_l = new TGraphErrors(1, x_l, y_l, ex_l, ey_l );
		
    gStyle->SetEndErrorSize(0);
    gr_l->SetMarkerSize(0.9);
    gr_l->Draw("0P");
		
    latex.SetTextFont(42);
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);    
    latex.SetTextSize(0.10);    
    latex.SetTextAlign(12); 
		
    TLine line_;
    TBox  box_;
    TBox  box_1;
    TBox  box_2;
    TBox  box_3;
    float xx_ = x_l[0];
    float yy_ = y_l[0];
    //Label title
    box_1.SetLineColor( histLineColor );
    box_1.SetFillColor( kYellow );
    box_1.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 );
    box_1.SetFillStyle(0);
    latex.DrawLatex(xx_+1.*bwx_,yy_,"EWK WW");
		
    yy_ -= gap_;
    box_.SetLineStyle( kSolid );
    box_.SetLineWidth( 1 );

    box_.SetLineColor( histLineColor );
    //Color of box
    box_.SetFillColor( kRed );
    box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 );
    box_.SetFillStyle(0);
    box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 );

    latex.DrawLatex(xx_+1.*bwx_,yy_,"EWK WZ");

    yy_ -= gap_;
    
    box_2.SetLineColor( histLineColor );
    box_2.SetFillColor( kBlue );
    box_2.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 );
    box_2.SetFillStyle(0);
    box_2.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 );

    latex.DrawLatex(xx_+1.*bwx_,yy_,"QCD WW");

    yy_ -= gap_;

    box_3.SetLineColor( histLineColor );
    //box_3.SetFillColor( );
    box_3.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 );
    box_3.SetFillStyle(0);
    box_3.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 );

    latex.DrawLatex(xx_+1.*bwx_,yy_,"W + jets");

    canv->cd();
  }

  {  

  //TFile file_("wwssjj_ewk.root","READ");
    //TFile *file_ewk = new TFile("wwssjj_ewk.root");
    //file_ewk->cd("VBSWW");
    
    //TH1D *h = (TH1D*)gDirectory->Get("mjj");
    //h->SetMarkerStyle(20);
    //h->SetMarkerSize(markerSize);

    //h->Draw("esamex0");

    //gDirectory->cd();

    TFile *file_wpj = new TFile("wwssjj_wpj_combined.root");
    file_wpj->cd("VBSWpJ");
    TH1D *h3 = (TH1D*)gDirectory->Get("mjj");
    h3->SetMarkerStyle(20);
    //h3->SetFillColor(1);
    h3->SetMarkerSize(markerSize);
    h3->Draw("histsame");
    
    gDirectory->cd();

    TFile *file_ewk = new TFile("wwssjj_ewk_1.root");
    file_ewk->cd("VBSWW_1");
    TH1D *h = (TH1D*)gDirectory->Get("mjj");
    h->SetMarkerStyle(20);
    h->SetFillColor(kYellow);
    h->SetMarkerSize(markerSize);
    h->Draw("histsame");

    gDirectory->cd();        

    TFile *file_qcd = new TFile("wwssjj_wzewk_1.root");
    file_qcd->cd("VBSWW_1");
    TH1D *h1 = (TH1D*)gDirectory->Get("mjj");
    h1->SetMarkerStyle(20);
    h1->SetFillColor(kRed);
    h1->SetMarkerSize(markerSize);
    h1->Draw("histsame");

    gDirectory->cd();

    TFile *file_wz = new TFile("wwssjj_qcd_1.root");
    file_wz->cd("VBSWW_1");
    TH1D *h2 = (TH1D*)gDirectory->Get("mjj");
    h2->SetMarkerStyle(20);
    h2->SetFillColor(kBlue);
    h2->SetMarkerSize(markerSize);
    h2->Draw("histsame");

  }

    CMS_lumi( canv, iPeriod, iPos );

    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();

    canv->Print(canvName+".pdf",".pdf");
    canv->Print(canvName+".png",".png");

    return canv;
}
