program CeliaWaterFlow;
//Coder : Ardiansyah (http://ardiansyah.net)
 
{$mode objfpc}{$H+}
 
uses
  Classes, SysUtils, Math
  { you can add units after this };
 
  const
  gr     = 9.8;      //gravitational constant
  im     = 0.000001; //maximum allowable mass balance error
  wd     = 1000;     //density of water (kg/m3)
  maxits = 200;      //maximum number of iterations
  mmax   = 100;      //maximum number of nodes
  R      = 8.3143;   //Gas constant (J/mol K)
  Mw     = 0.018;    //molecular mass of water (kg/mol)
  //Dv     = 0.000024; //vapor diffusivity (m2/s)
  eps    = 0.5;
 
type
  onedims = array [0..mmax] of single;
  onedimd = array [0..mmax] of double;
 
var
  p   : onedimd;                                                                //node water potential at start of time step
  pi  : onedimd;
  pn  : onedimd;
  z   : onedims;                                                                //node depths
  cpu : onedims;                                                                //upper node water capacity
  cpl : onedims;                                                                //lower node water capacity
  v   : onedims;                                                                //soil volume at node
  wu  : onedims;                                                                //upper water content in middle of time step
  wl  : onedims;                                                                //lower water content in middle of time step
  wnu : onedims;                                                                //upper water content at end of time step
  wnl : onedims;                                                                //lower water content at end of time step
  wiu : onedims;                                                                //upper water content at initial timestep
  wil : onedims;						  										//lower water content at initial timestep
  w   : onedims;                                                                //water content to be plot
  Jl  : onedimd;                                                                //element liquid fluxes
  Jv  : onedimd;                                                                //element vapor flux, for further coding
  ks  : onedims;                                                                //element saturated conductivity (cm/s)
  ws  : onedims;                                                                //element saturated water content
  pe  : onedimd;                                                                //element air entry potential (mm H2O)
  dwudp,dwldp,dwnudp,dwnldp :onedims; //derrivative of w over p
  //h   : onedimd;                                                              //node humidity
 
  m   : integer;                         //node
  //PET,ha       : single;                                                      //potential evaporation rate, atmospheric humidity
  suminf       : single;                                                        //cumulative infiltration if there is exsist infiltration
  File2        : TextFile;
  depth, bd, wp, ksat, ae, wsat, time, dt, endtime : single;
 
function Power(Base, Exponent: single): single;
begin
  Power := Exp(Exponent*ln(Base));
end;{end of power}
 
procedure GeoGrid(depth:single; m:integer);
var sc, dz :single;
    i      :integer;
begin
  sc := 0;
  for i := 1 to m do sc := sc + i*i;
  dz := depth/sc;
  z[0] := 0;
  z[1] := 0;
  for i := 1 to m do z[i+1] := z[i] + dz*i*i;
end;
 
procedure LinGrid(depth:single; m:integer);
var i : integer;
    dz: single;
begin
  dz   := depth/m;
  z[0] := 0;
  z[1] := 0;
  for i:= 1 to m do z[i+1] := z[i]+dz;
end;
 
procedure Average_w(m: integer; var w : onedims);
var i : integer;
begin
  for i := 1 to m+1 do w[i] := 0.5 * (wnu[i] + wnl[i]);
end;
 
function WaterCapacity(p:single):single;
var a,b,c,d,e,f : single;
begin
  a := 0.343;
  b := 0.04133;
  c := 5.603;
  d := 0.8215;
  e := 0.08624;
  f := 14.37;
  p := abs(p);
  WaterCapacity := abs( (-a*d/power((1+power((b*p),c)),(d+1)))
                   *((c*power((b*p),c))/p) - e/(f*(p+1)) );
end;
 
procedure WaterCont(p :single; var w, dwdp : single);
var a,b,c,d,e,f : single;
begin
  a := 0.343;
  b := 0.04133;
  c := 5.603;
  d := 0.8215;
  e := 0.08624;
  f := 14.37;
  p := abs(p);
  w := a/power((1+power((b * p),c)),d)  +  e*(1-(Log10(p+1)/f));
  dwdp := WaterCapacity(p);
end;{WaterCont}
 
procedure HydrConductivity(ks,wc,ws :single; var kh :single);
//hydraulic conductivity is for layer or between node
var a1,b1,c1,x1,x2 : single;
begin
  a1 := 40.3;
  b1 := 0.672;
  c1 := 3.31;
  x1    := power((wc/ws),b1);
  x2    := 1-x1;
  If (wc = ws) then kh := ks else kh := ks * Exp(-a1*(power(x2,c1)));
end;{HydrConductivity}
 
procedure k_bar(i : integer;var ki, kn :single);
//k average in time, ingat water content di node bukan di layer
//wat content di node, tapi di layer disubscript sebagai i+1/2
var wi, wn, kbar : single;
const eps = 0.5;
begin
  wi   := 0.5*(wiu[i]+wil[i]); //upper at i itu = i, lower at i itu = i+1
  wn   := 0.5*(wnu[i]+wnl[i]);
  HydrConductivity(ks[i],wi,ws[i],ki);
  HydrConductivity(ks[i],wn,ws[i],kn);
  kbar := (1-eps)*ki + eps*kn; //time averaged k, not used in this calculation
end;
 
procedure InitSoil(m:integer; wp, ae, ksat, wsat : single);
var i : integer;
begin
  for i := 1 to m + 1 do
  begin
    ks[i]:= ksat;                                                               // p.56, Campbell
    ws[i]:= wsat;
    pe[i]:= ae;                                                                 // Air entry potential
    p[i] := wp;                                                                 // Water potential
    WaterCont(p[i],wu[i],dwudp[i]);
    wl[i]:= wu[i];
    v[i] := (z[i+1]-z[i-1])/2; //dz
  end;
  z[0]    := -1E+20;//0;//
  z[m+1]  := +1E+20;                                                            // No upward vapor flux into bottom
  wiu     := wu;
  wil     := wl;
  wnu     := wu;
  wnl     := wl;
  wnl[m+1]:= wl[m];
  wnu[m+1]:= wu[m];
  dwnudp  := dwudp;
  p[0]    := p[1];
  pi      := p;     //pi initial timestep, p middle of timestep, pn akhir timestep
  pn      := p;
  suminf  := 0;
end;
 
procedure q_liquid(i:integer; var qi, qn, J_liquid :double);
var ki,kn : single;
//  'i is element number; upper p is i, lower is i+1
begin
  k_bar(i,ki,kn);
  qi := -(ki/(z[i+1]-z[i])*(pi[i+1]-pi[i]))+ki;
  qn := -(kn/(z[i+1]-z[i])*(pn[i+1]-pn[i]))+kn;
  J_liquid := (1-eps)*qi + eps*qn;                  //flux hasil perhitungan antara dua timestep, bukan dua level iterasi
end;
 
procedure LiqCoeff(i:integer; var UpCoefn,LowCoefn,ResCoefn:single);
var wn, kn : single;
begin
  wn := 0.5*(wnu[i]+wnl[i]);
  HydrConductivity(ks[i],wn,ws[i],kn);
  UpCoefn  := eps*(kn/(z[i+1]-z[i]));
  LowCoefn := -eps*(kn/(z[i+1]-z[i]));
  ResCoefn := eps*kn;
end;
 
procedure SolveEvap(m :integer; dt, flux, psurface :single; success :boolean; nits :integer);
var  i                           : integer;
     se, sw                      : single;
     A,Bx,C,D,qli,qln            : onedimd;
     UpCoefn, LowCoefn, ResCoefn : onedims;		//coefficient for a,b,c and d
 
begin
  Nits := 0;
  //Apply upper boundary condition
  If psurface &lt; 0 then   begin     If psurface &gt; pe[1] then p[1] := pe[1] else p[1] := psurface;
  end
  else
  begin
    Jl[0]  := -flux; {flux boundary condition}
    qli[0] := Jl[0]; //time average flux for flux boundary condition
    qln[0] := Jl[0];
  end;
 
  Repeat
    se   := 0;                                                                  // Sum of error
    Nits := Nits + 1;                                                           // Number of iteration
    p := pn; // untuk evaluasi error
    for i := 1 to m do  //jika i mulai dari layer kedua karena i=1 untuk boundary condition
    begin
      cpu[i] := v[i]*dwnudp[i]/dt;  // Water capacity
      q_liquid(i, qli[i],qln[i],Jl[i]);
      LiqCoeff(i,UpCoefn[i],LowCoefn[i],ResCoefn[i]);
 
      A[i] := UpCoefn[i-1];
      C[i] := -LowCoefn[i];
      Bx[i]:= -UpCoefn[i]+LowCoefn[i-1] - cpu[i];
      if i = 1 then
      begin
        D[i] := -(1-eps)*(qli[i]-qli[i-1]) - cpu[i]*p[i]+(wu[i]-wiu[i])*(v[i]/dt) + ResCoefn[i]-ResCoefn[i-1]+eps*qln[i-1]; //involve initial state and updated iteraton level
      end
      else
      begin
        D[i] := -(1-eps)*(qli[i]-qli[i-1]) - cpu[i]*p[i]+(wu[i]-wiu[i])*(v[i]/dt) + ResCoefn[i]-ResCoefn[i-1]; //involve initial state and updated iteraton level
      end;
    end;
 
    //Thomas Algorithm starts here
    If psurface &lt; 0 then //dirichlet BC not flux BC, start calculation from i = 2
    begin
      D[1] := 0;
      C[1] := 0;
    end;{if}
    for i := 1 to m - 1 do
    begin
      C[i] := C[i]/Bx[i];
      D[i] := D[i]/Bx[i];
      Bx[i+1] := Bx[i+1] - A[i+1]*C[i];
      D[i+1]  := D[i+1] - A[i+1]*D[i];
    end;{for}
    pn[m] := D[m]/Bx[m];
    for i := m - 1 downto 1 do
    begin
      pn[i] := D[i] - C[i] * pn[i+1];
    end;{for}
    //end thomas algorithm
 
    //calculate water content at each node and total error of mass balance
    for i := 1 to m do
    begin
      WaterCont(p[i],wu[i],dwudp[i]);
      WaterCont(p[i+1],wl[i],dwldp[i]);
      WaterCont(pn[i],wnu[i],dwnudp[i]);
      WaterCont(pn[i+1], wnl[i], dwnldp[i]);
      se := se + abs(wu[i]-wnu[i]);
    end;
  Until (se &lt; im) Or (Nits &gt; maxits);
 
  pn[m+1] := pn[m];                                                               // Unit gradient drainage at bottom, no drainage
  If nits &lt;= maxits then   begin     sw := 0;  //water storage     for i := m downto 1 do     begin       sw    := sw + v[i] * (wnu[i] - wiu[i] + wnl[i-1] - wil[i-1])/2;//cm3 water /cm2 soil     end;{for}     wu  := wnu;     wl  := wnl;     wiu := wnu;     wil := wnl;     pi  := pn;   //memperkenalkan "pi" untuk assign pn ke p timestep sebelumnya yang bisa digunakan dalam menghitung flux     //flux := Jv[0];//Jv sebagai flux digunakan untuk kasus evaporasi dan evaporasi berubah2 tiap timestep, bukan infiltrasi, kalo kasus infiltrasi, flux = Jl[0]     suminf := suminf + flux*dt;     success:= True;   end   else success := False;   Average_w(m,w); end; var  i, nits : integer;      pub     : single;                                                          //water potential at the upper boundary; if pub&gt;0 then upper bc is flux
     flux    : single;                                                          //flux of water into soil (cm/s)
     success : boolean;
     //variable di unit 2
     day     : integer;
{$IFDEF WINDOWS}{$R CeliaWaterFlow.rc}{$ENDIF}
 
begin
//Soil Parameter
  m       := 15;                                             // Number of elements
  depth   := 50;                                             // Depth (cm)
  bd      := 1.4;                                            // Soil bulk density (Mg/m3 = g/cm3)
  wp      := -10.2;                                          // Initial water potential (cmH2O)
  ae      := -10.2;                                          // Air entry water potential (cmH2O)
  ksat    := 0.01093;                                        // Saturated conductivity (cm/s)
  wsat    := 0.43;                                           // Saturated water content (cm3/cm3)
  dt      := 600;                                            // Time step (s)
  endtime := 24;                                             // Length of simulation (hr)
  time    := 0;                                              // Start time
 
  If ae &gt; 0 then ae := -ae;
  LinGrid(depth,m);//GeoGrid(depth, m);//                                                        Discretization, whether grid is linear or geometric
  day  := 0;
 
  InitSoil(m, wp, ae, ksat, wsat);
  pub  := 1;                                                 // p at surface =1, Set upper boundary condition for evaporation (positive)
                                                             // pub &gt;=0 flux boundary condition
 
  //create text file
  Assign(File2,'/home/ardiansyah/Desktop/CeliaWatFlow.txt');
  Rewrite(File2);
  Write(File2,' node '); for i:= 1 to m do Write(File2,i,' '); Writeln(File2,'');
  Write(File2,' Depth(m) '); for i:= 1 to m do Write(File2,z[i],' '); Writeln(File2,'');
  Writeln(File2,' DOY ', ' time',' Wat.Cont(m3/m3)');
 
  Repeat
    time := time + dt/3600;
    flux := -0.000009259; //0.000009259 cm/s = 8 mm/day, load flux into surface  (-) means upward flow
    SolveEvap(m,dt,flux,pub,success,nits);
 
    Write(File2,{daystrt,} time);
    for i := 1 to m do
    begin
      Write(File2,wnu[i],' ');                                                    // Print to file
      Sleep(5);
    end;{for}
    Write(time);
    Writeln(wnu[1]);
    Writeln(File2,' ');
 
  Until time&gt;endtime;
  Close(File2);
  Readln;
end.
