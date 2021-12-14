#include <iostream>
#include <stdlib.h>
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TFile.h"

void Format(TH1F* Hist, int LineWidth, int LineColor, int FillStyle, int FillColor)
{
	Hist->SetLineWidth(LineWidth);
	Hist->SetLineColor(LineColor);
	Hist->SetFillStyle(FillStyle);
	Hist->SetFillColor(FillColor);
}





void Titleold(TH1F* Hist, const char *Title, const char *XTitle, const char *YTitle, float TitleOffset)
{
	Hist->SetTitle(Title);
	Hist->GetXaxis()->SetTitle(XTitle);
	Hist->GetXaxis()->CenterTitle();	
	// Hist->GetXaxis()->SetTitleOffset(TitleOffset);
	Hist->GetYaxis()->SetTitle(YTitle);
	Hist->GetYaxis()->CenterTitle();
	Hist->GetYaxis()->SetTitleOffset(TitleOffset);
}

void Title(TH1F* Hist, const char *Title, const char *XTitle, const char *YTitle, float TitleOffset, float TitleSize, float LabelSize)
{
    Hist->SetTitle(Title);
    Hist->GetXaxis()->SetTitle(XTitle);
    //Hist->GetXaxis()->SetTitleSize(0.055);
    Hist->GetXaxis()->SetTitleSize(TitleSize);
    //Hist->GetXaxis()->CenterTitle();  
    Hist->GetXaxis()->SetTitleOffset(TitleOffset);
    Hist->GetXaxis()->SetLabelSize(LabelSize);
    Hist->GetYaxis()->SetTitle(YTitle);
    //Hist->GetYaxis()->CenterTitle();
    Hist->GetYaxis()->SetTitleOffset(TitleOffset+0.1);
    Hist->GetYaxis()->SetTitleSize(TitleSize);
    Hist->GetYaxis()->SetLabelSize(LabelSize);
}




void histo4_cut2(std::string outfilename) {
    
 gStyle->SetCanvasPreferGL(kTRUE);


   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",1562,75,1000,1000);
   //  TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",1562,75,1154,698);   //old dimension

   Canvas_1->Range(-12.5,-1.25,112.5,11.25);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetTicks(1);

   //Open the root file
TFile *f1 = TFile::Open("finalcuts3.root");	
TDirectory* directory = (TDirectory*)f1->Get("analyzeBasicPat");
//TH1F *DEN1, *TRIG0, *TRIG1, *TRIG2, *DEN2, *TRIG02, *TRIG12, *TRIG22 ;
TH1F *histo1, *histo2, *histo3, *histo4; 

//auto c2h = new TCanvas("c2h","c2h",10,10,800,600);

directory->GetObject("metPt_case2", histo1);
directory->GetObject("metPt_case2_trigger0",histo2);
directory->GetObject("metPt_case2_trigger1",histo3);
directory->GetObject("metPt_case2_trigger2",histo4);

/*
directory->GetObject("metPt_case1", DEN2);
directory->GetObject("metPt_case1_trigger0",TRIG02);
directory->GetObject("metPt_case1_trigger1",TRIG12);
directory->GetObject("metPt_case1_trigger2",TRIG22);
*/

Canvas_1->cd();
histo1->Rebin(4);
histo2->Rebin(4);
histo3->Rebin(4);
histo4->Rebin(4);
histo1->SetName("Cuts 2");
histo2->SetName("Trigger 1 + Cuts 2");
histo3->SetName("Trigger 2 + Cuts 2");
histo4->SetName("Trigger 3 + Cuts 2");
//NUM->Rebin(4);
//DEN->Draw();
//DEN->SetLineColor(kRed);
//DEN->Draw();
//NUM->Draw("sames");

TH1F *RATIO1 = (TH1F*)histo2->Clone("RATIO1");
RATIO1->Divide(histo1); 

TH1F *RATIO2 = (TH1F*)histo3->Clone("RATIO2");
RATIO2->Divide(histo1); 

TH1F *RATIO3 = (TH1F*)histo4->Clone("RATIO3");
RATIO3->Divide(histo1); 

//auto gr = new TGraphAsymmErrors();
//  gr->Divide(histo2,histo1,"cl=0.683 b(1,1) mode");
  auto gr1 = new TGraphAsymmErrors(histo2,histo1,"cl=0.683 b(1,1) mode");
    auto gr2 = new TGraphAsymmErrors(histo3,histo1,"cl=0.683 b(1,1) mode");
      auto gr3 = new TGraphAsymmErrors(histo4,histo1,"cl=0.683 b(1,1) mode");
        
        Title(histo1, " ","Title X" ,"Title Y",  0.85, 0.055, 0.05);
        Title(histo2, " ","MET [GeV]" ,"Entries", 0.85, 0.055, 0.03);
        Title(histo3, " ","Title X" ,"Title Y", 0.85, 0.055, 0.05);
        Title(histo4, " ","Title X" ,"Title Y", 0.85, 0.055, 0.05);

  
      Format(histo2, 2, 634, 1001, 623);    
      Format(histo1, 2, 1, 0, 591);
  
    
    
      Format(histo3, 2, 420,1001, 419);
      Format(histo4, 2, 602, 1001, 591);   

        histo2->SetFillColorAlpha(623, 0.5);
        histo3->SetFillColorAlpha(419, 0.5);
        histo4->SetFillColorAlpha(591, 0.5);

        
    // StatBox
       // histo2->GetXaxis()->SetRangeUser(0, 400);
        histo2->Draw();
    	gPad->Update();
		TPaveStats* St2 = (TPaveStats*) histo2->FindObject("stats");
    St2->SetX1NDC(0.52);
		St2->SetX2NDC(0.68);
		St2->SetY1NDC(0.75);
		St2->SetY2NDC(0.88);
    
        St2->SetBorderSize(1);
        St2->SetFillColor(0);
        St2->SetTextAlign(12);
        St2->SetTextColor(634);
        St2->SetTextFont(42);
    

    
    // StatBox
     histo1->Draw("sames");
     histo3->Draw("sames");
     histo4->Draw("sames");
    
    // Defining the proper maximum range for the histogram
       
        float maxhisto1=histo1->GetMaximum();
        float maxhisto2=histo2->GetMaximum();
        float maxhisto3=histo3->GetMaximum();
        float maxhisto4=histo4->GetMaximum();
    float max1;
    float max2;
    
        if(maxhisto1 > maxhisto2){
        max1=maxhisto1;
                                }
        else
            {
            max1=maxhisto2;
            }
    
        if(maxhisto3 > maxhisto4){
        max2=maxhisto3;
                                }
        else
            {
            max2=maxhisto4;
            }
    
          if(max1 > max2){
         histo2->SetMaximum(max1*1.2);
                                }
        else
            {
         histo2->SetMaximum(max2*1.2);
            }
    
    
        		gPad->Update();
		TPaveStats* St1 = (TPaveStats*) histo1->FindObject("stats");
    St1->SetX1NDC(0.34);
		St1->SetX2NDC(0.5);
		St1->SetY1NDC(0.75);
		St1->SetY2NDC(0.88);
       
        St1->SetBorderSize(1);
        St1->SetFillColor(0);
        St1->SetTextAlign(12);
        St1->SetTextColor(1);
        St1->SetTextFont(42);
  
  
    
    
    		gPad->Update();
		TPaveStats* St3 = (TPaveStats*) histo3->FindObject("stats");


    St3->SetX1NDC(0.34);
		St3->SetX2NDC(0.5);
		St3->SetY1NDC(0.60);
		St3->SetY2NDC(0.73);
 
       
        St3->SetBorderSize(1);
        St3->SetFillColor(0);
        St3->SetTextAlign(12);
        St3->SetTextColor(420);
        St3->SetTextFont(42);
    
    		gPad->Update();
		TPaveStats* St4 = (TPaveStats*) histo4->FindObject("stats");
		St4->SetX1NDC(0.52);
		St4->SetX2NDC(0.68);
		St4->SetY1NDC(0.60);
		St4->SetY2NDC(0.73);
       
        St4->SetBorderSize(1);
        St4->SetFillColor(0);
        St4->SetTextAlign(12);
        St4->SetTextColor(602);
        St4->SetTextFont(42);


/*
  TPaveText *pt4 = new TPaveText(0.0681363   ,0.872464    ,0.691383  ,0.988406  ,"blNDC");
    pt4->AddText("General title Histo");
    pt4->SetBorderSize(0);
    pt4->SetFillColor(0);
    pt4->SetFillStyle(0);
    pt4->SetTextFont(42);
    pt4->SetTextAlign(11);
    pt4->SetTextSize(0.08);
    pt4->Draw();


     TPaveText *pt5 = new TPaveText(0.131263  ,0.655072    ,0.232465   ,0.992754   ,"blNDC");
    pt5->AddText("CMS");
   // pt1->AddText("1.63 ");
   // pt1->AddText("1.67 ");
    pt5->SetBorderSize(0);
    pt5->SetFillColor(0);
    pt5->SetFillStyle(0);
    pt5->SetTextFont(62);
    pt5->SetTextAlign(11);
    //pt5->GetTextSize();
    pt5->Draw();


        TPaveText *pt6 = new TPaveText(0.103206   ,0.746377   ,0.726453   ,0.805797    ,"blNDC");
    pt6->AddText("Preliminary");
   // pt1->AddText("1.63 ");
   // pt1->AddText("1.67 ");
    pt6->SetBorderSize(0);
    pt6->SetFillColor(0);
    pt6->SetFillStyle(0);
    pt6->SetTextFont(52);
    pt6->SetTextAlign(11);
    //pt6->GetTextSize();
    pt6->Draw();

    */

        TLegend *legend = new TLegend(0.726453 ,0.697128,0.875752  ,0.856397,NULL,"brNDC");
        legend->SetBorderSize(0);
        legend->AddEntry(histo1, "Cuts 2","f");
        legend->AddEntry(histo2, "Trigger 1 + Cuts 2","f");
        legend->AddEntry(histo3, "Trigger 2 + Cuts 2","f");
        legend->AddEntry(histo4, "Trigger 3 + Cuts 2","f");
        legend->Draw();
    



    
    TFile *out = new TFile(outfilename.c_str(), "RECREATE");
    

if (out->IsZombie()) { cout<<"Open not succesful"<<endl;}
else { cout<<outfilename.c_str()<<" succesfully created"<<endl;
  histo2->Write();
    histo1->Write();
    Canvas_1->Write();
    histo3->Write();
    histo4->Write();
    RATIO1->Write();
RATIO2->Write();
RATIO3->Write();

gr1->Write();
gr2->Write();
gr3->Write();


}

 
    
  
    
     Canvas_1->Update();
	 Canvas_1->SaveAs("histo4cut2a.png");
    
    

}
