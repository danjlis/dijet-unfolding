
std::vector<float> dividebin(TH2* h2, int bin, int ndivs)
{
  std::vector<float> divisions;
  int binx, biny, binz;
  h2->GetBinXYZ(bin, binx, biny, binz);
  float x1 = h2->GetXaxis()->GetBinLowEdge(binx);
  float y1 = h2->GetYaxis()->GetBinLowEdge(biny);
  float x2 = x1 +  h2->GetXaxis()->GetBinWidth(binx);
  float y2 = y1 +  h2->GetYaxis()->GetBinWidth(biny);
  float total_area = (x2 - x1) * (y2 - y1);
  // std::cout << x1 << " , " << y1 << " -->  " << x2 << " , " << y2 << " = " << total_area << std::endl;
  float area_total = 0;
  float dx = 1./(float)ndivs;
  // std::cout << dx << std::endl;
  for (int i = 0; i <= ndivs; i++)
    {
      float slope = ((float)i)*dx;
      // std::cout << "slope =  " << slope << std::endl;
      float yb1 = slope*x1;
      float yb2 = slope*x2;
      // std::cout << "YBin : " << yb1 << " -- " << yb2 << std::endl;
      if (yb1 <= y1 && yb2 <= y1)
	{
	  divisions.push_back(0);
	  continue;
	}
      if (yb1 >= y2 && yb2 > y2)
	{
	  divisions.push_back((total_area - area_total)/total_area);
	  area_total = 25;
	  continue;
	}

      float xb1 = y1/slope;
      float xb2 = y2/slope;

      // std::cout << "XBin : " << xb1 << " -- " << xb2 << std::endl;

      if (xb1 < x1) xb1 = x1;
      if (xb2 > x2) xb2 = x2;

      if (yb1 < y1) yb1 = y1;
      if (yb2 > y2) yb2 = y2;
      // std::cout << "YBin : " << yb1 << " -- " << yb2 << std::endl;
      // std::cout << "XBin : " << xb1 << " -- " << xb2 << std::endl;
      

      float area = 0.5*(xb2 - xb1)*(yb2 - yb1);
      if (yb1 > y1)
	{
	  area += (yb1 - y1)*(x2 - x1);
	}
      if (xb2 < x2)
	{
	  area += (y2 - yb1)*(x2 - xb2);
	}

      area -= area_total;
      area_total += area;
      
      // std::cout << area << std::endl;
      divisions.push_back(area/total_area);
    }
  //  divisions.push_back(total_area - area_total);
  return divisions;  
}

void dividebins()

{
  TH2D *h2=new TH2D("h2","", 10, 0, 50, 10, 0, 50);

  for (int i = 1; i <= h2->GetNbinsX(); i++)
    {
      for (int j = 1; j <= i ; j++)
	{
	  int bin = h2->GetBin(i,j);
	  std::vector<float> divs = dividebin(h2, bin, 10);
	  std::cout << i << " / " << j << ": ";
	  float total = 0;
	  for (int id = 0; id < divs.size();id++)
	    {
	      total += divs.at(id);
	    std:cout << divs.at(id) << " + ";
	    }
	  std::cout << " = " << total<< std::endl;
	}

    }
}
