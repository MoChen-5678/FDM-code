!*****************************************************************************!
!
Module jsymbols
!
!*****************************************************************************!

IMPLICIT NONE

    integer, parameter, private :: igfv = 100
    integer, dimension(0:igfv), private :: iv
    double precision, dimension(0:igfv), private :: sq, sqi, sqh, shi, gmi, wg
    double precision, dimension(0:igfv), private :: fak, fad, fi, wf, wfi, gm2
    double precision, dimension(0:igfv), private :: wgi

!--- Define math and physics constants    
    double precision, private :: one, two, half, third, zero, pi
    
    data one/1.0d0/, two/2.0d0/, half/0.5d0/, zero/0.0d0/
    data third/0.3333333333333333333333333333333333333333333333333333333333d0/
    data pi/3.141592653589793d0/
!
!--- This module contains the subroutines which calculate the 
!--- 3j symbol, 6j symbol and 9j symbol
!
CONTAINS
!=============================================================================!
!
SUBROUTINE neufj(xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33,c9j)
!
!=============================================================================!
!
!--- Subroutine for 9j symbol
!---
    double precision, intent(in)    :: xj11, xj12, xj13, xj21, xj22, xj23
    double precision, intent(in)    :: xj31, xj32, xj33
    double precision, intent(out) :: c9j

    double precision, dimension(301) ::flog
    double precision :: eps1, eps2
    integer :: i
    data(flog(i),i=2,31)/ 0.d0,.69314718d0,1.7917595d0,3.1780538d0, &
     & 4.7874917d0,6.5792511d0,8.5251613d0,10.604603d0,12.801827d0, &
     & 15.104413d0,17.502307d0,19.987214d0,22.552163d0,25.191221d0, &
     & 27.899271d0,30.671860d0,33.505072d0,36.395445d0,39.339884d0, &
     & 42.335616d0,45.380139d0,48.471180d0,51.606674d0,54.784729d0, &
     & 58.003604d0,61.261702d0,64.557537d0,67.889743d0,71.257038d0, &
     & 74.658235d0/
    data(flog(i),i=32,61)/78.092223d0,81.557959d0,85.054466d0,88.580827d0, &
     & 92.136175d0,95.719694d0,99.330612d0,102.96820d0,106.63176d0, &
     & 110.32064d0,114.03421d0,117.77188d0,121.53308d0,125.31727d0, &
     & 129.12393d0,132.95257d0,136.80272d0,140.67392d0,144.56574d0, &
     & 148.47776d0,152.40959d0,156.36083d0,160.33112d0,164.32011d0, &
     & 168.32744d0,172.35279d0,176.39584d0,180.45629d0,184.53383d0, &
     & 188.62817d0/
    data(flog(i),i=62,91)/192.73904d0,196.86618d0,201.00931d0,205.16820d0, &
     & 209.34258d0,213.53224d0,217.73693d0,221.95644d0,226.19054d0, &
     & 230.43904d0,234.70172d0,238.97839d0,243.26885d0,247.57291d0, &
     & 251.89040d0,256.22113d0,260.56494d0,264.92164d0,269.29110d0, &
     & 273.67312d0,278.06757d0,282.47429d0,286.89313d0,291.32394d0, &
     & 295.76659d0,300.22094d0,304.68685d0,309.16419d0,313.65283d0, &
     & 318.15264d0/
    data(flog(i),i=92,121)/322.66349d0,327.18529d0,331.71788d0,336.26118d0,&
     & 340.81505d0,345.37940d0,349.95411d0,354.53908d0,359.13420d0, &
     & 363.73937d0,368.35449d0,372.97946d0,377.61419d0,382.25859d0, &
     & 386.91255d0,391.57598d0,396.24881d0,400.93094d0,405.62230d0, &
     & 410.32277d0,415.03230d0,419.75080d0,424.47819d0,429.21439d0, &
     & 433.95932d0,438.71291d0,443.47508d0,448.24576d0,453.02489d0, &
     & 457.81238d0/
    data(flog(i),i=122,151)/462.60817d0,467.41220d0,472.22438d0,     &
     & 477.04466d0,481.87298d0,486.70926d0,491.55345d0,496.40547d0, &
     & 501.26529d0,506.13282d0,511.00802d0,515.89082d0,520.78117d0, &
     & 525.67901d0,530.58428d0,535.49694d0,540.41692d0,545.34417d0, &
     & 550.27865d0,555.22029d0,560.16905d0,565.12488d0,570.08772d0, &
     & 575.05753d0,580.03427d0,585.01787d0,590.00830d0,595.00552d0, &
     & 600.00946d0,605.02010d0/
    data(flog(i),i=152,181)/610.03738d0,615.06126d0,620.09170d0,     &
     & 625.12866d0,630.17208d0,635.22193d0,640.27818d0,645.34077d0, &
     & 650.40968d0,655.48486d0,660.56626d0,665.65385d0,670.74760d0, &
     & 675.84747d0,680.95341d0,686.06541d0,691.18340d0,696.30735d0, &
     & 701.43726d0,706.57306d0,711.71472d0,716.86221d0,722.01551d0, &
     & 727.17456d0,732.33934d0,737.50983d0,742.68598d0,747.86776d0, &
     & 753.05516d0,758.24811d0/
    data(flog(i),i=182,211)/763.44661d0,768.65061d0,773.86010d0,     &
     & 779.07503d0,784.29539d0,789.52114d0,794.75224d0,799.98869d0, & 
     & 805.23044d0,810.47747d0,815.72973d0,820.98722d0,826.24991d0, &
     & 831.51778d0,836.79078d0,842.06890d0,847.35209d0,852.64036d0, &
     & 857.93366d0,863.23199d0,868.53529d0,873.84356d0,879.15676d0, &
     & 884.47488d0,889.79789d0,895.12577d0,900.45848d0,905.79603d0, &
     & 911.13836d0,916.48547d0/
    data(flog(i),i=212,241)/921.83732d0,927.19391d0,932.55521d0,     &
     & 937.92118d0,943.29181d0,948.66710d0,954.04699d0,959.43148d0, &
     & 964.82056d0,970.21419d0,975.61235d0,981.01503d0,986.42220d0, &
     & 991.83385d0,997.24995d0,1002.6705d0,1008.0954d0,1013.5248d0, &
     & 1018.9585d0,1024.3966d0,1029.8389d0,1035.2857d0,1040.7367d0, &
     & 1046.1920d0,1051.6516d0,1057.1155d0,1062.5836d0,1068.0558d0, &
     & 1073.5323d0,1079.0129d0/     
    data(flog(i),i=242,271)/1084.4977d0,1089.9866d0,1095.4797d0,     &
     & 1100.9768d0,1106.4781d0,1111.9834d0,1117.4928d0,1123.0063d0, &
     & 1128.5237d0,1134.0452d0,1139.5706d0,1145.1001d0,1150.6335d0, &
     & 1156.1708d0,1161.7120d0,1167.2573d0,1172.8063d0,1178.3593d0, &
     & 1183.9161d0,1189.4768d0,1195.0413d0,1200.6097d0,1206.1818d0, &
     & 1211.7577d0,1217.3375d0,1222.9209d0,1228.5082d0,1234.0992d0, &
     & 1239.6939d0,1245.2924d0/
    data(flog(i),i=272,301)/1250.8944d0,1256.5003d0,1262.1097d0,     &
     & 1267.7228d0,1273.3396d0,1278.9600d0,1284.5840d0,1290.2117d0, &
     & 1295.8429d0,1301.4777d0,1307.1160d0,1312.7580d0,1318.4034d0, &
     & 1324.0524d0,1329.7048d0,1335.3609d0,1341.0203d0,1346.6833d0, &
     & 1352.3497d0,1358.0196d0,1363.6929d0,1369.3697d0,1375.0499d0, &
     & 1380.7334d0,1386.4204d0,1392.1107d0,1397.8045d0,1403.5016d0, &
     & 1409.2020d0,1414.9058d0/
    data eps1,eps2/.1d0,-.2d0/
    
    double precision :: xn, xj, xjmin, xjmax, f, s, f1, s1
    double precision :: f2, s2, f3, s3, s2n, s2d, c2n, c2d
    
    integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12
    integer :: n13, n14, n15, n16, n17, n18, n19, n20, n21, n22
    integer :: n23, n24, n25, n26, n27, n28, n29, n30, n31, n32
    integer :: n33, n34, n35, n36, n37, n38, n39, n40, n41, n42
    integer :: n43, n44, n45, n46, n47, n48, n49, n50, n51, n52
    integer :: n53, n54, n55, n56, n57, n58, n59
    
    integer :: k, n, l, k1, l1, i1, i1m1, nn1, nd1, l2, i2, k2
    integer :: i2m2, nn2, nd2, k3, l3, i3, i3m3, nn3, nd3, km1
    integer :: kp1, kp2, kp3, km2, km3, jmin, jmax


!--- Calcul des combinaisons xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33
!                                                                                               
    xn   = -xj11 + xj21 + xj31 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n1   =  xn
    xn   = -xj32 + xj33 + xj31 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n2   =  xn
    xn   =  xj11 - xj21 + xj31 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n5   =  xn
    xn   =  xj32 - xj33 + xj31 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n6   =  xn
    xn   =  xj11 + xj21 - xj31 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n9   =  xn
    xn   =  xj32 + xj33 - xj31 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n10  =  xn
    xn   =  xj11 + xj32 + xj21 + xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n15  = xn
    n15  = n15 + 1
    xn   = xj11 + xj21 + xj31 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n16  = xn
    n16  = n16 + 1
    xn   = xj32 + xj33 + xj31 + eps1   
    if(xn.lt.0.d0) xn  = xn + eps2
    n17  = xn    
    n17  = n17 + 1  
    xn   =  - xj12 + xj22 + xj32 + eps1  
    if(xn.lt.0.d0) xn  = xn + eps2
    n21  = xn    
    xn   = - xj21 + xj22 + xj23 + eps1  
    if(xn.lt.0.d0) xn  = xn + eps2
    n23  = xn    
    xn   = xj12 - xj22 + xj32 + eps1   
    if(xn.lt.0.d0) xn  = xn + eps2
    n25  = xn    
    xn   = xj21 - xj22 + xj23 + eps1   
    if(xn.lt.0.d0) xn  = xn + eps2
    n27  = xn    
    xn   = xj12 + xj22 - xj32 + eps1   
    if(xn.lt.0.d0) xn  = xn + eps2
    n29  = xn    
    xn   = xj21 + xj22 - xj23 + eps1   
    if(xn.lt.0.d0) xn  = xn + eps2
    n31  = xn
    xn   =  - xj12 - xj21 + xj32 + xj23 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n33  = xn
    xn   = xj12 + xj22 + xj32 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n36  = xn
    n36  = n36 + 1
    xn   = xj21 + xj22 + xj23 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n38  = xn
    n38  = n38 + 1
    xn   =  - xj13 + xj23 + xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n41  = xn
    xn   =  - xj13 + xj11 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n44  = xn
    xn   = xj13 - xj23 + xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n45  = xn
    xn   = xj13 - xj11 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n48  = xn
    xn   = xj13 + xj23 - xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n49  = xn
    xn   = xj13 + xj11 - xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n52  = xn
    xn   =  - xj23 - xj11 + xj33 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n54  = xn
    xn   = xj13 + xj23 + xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n56  = xn
    n56  = n56 + 1
    xn   = xj13 + xj11 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n59  = xn
    n59  = n59 + 1
!
!--- test sur les xj11,xj12,xj13,xj21,xj22,xj23,xj31,xj32,xj33
!
    if(n9.lt.0) go to 50
    if(n5.lt.0) go to 50
    if(n1.lt.0) go to 50
    if(n10.lt.0) go to 50
    if(n6.lt.0) go to 50
    if(n2.lt.0) go to 50
    if(n29.lt.0) go to 50
    if(n25.lt.0) go to 50
    if(n21.lt.0) go to 50
    if(n31.lt.0) go to 50
    if(n27.lt.0) go to 50
    if(n23.lt.0) go to 50
    if(n49.lt.0) go to 50
    if(n45.lt.0) go to 50
    if(n41.lt.0) go to 50
    if(n52.lt.0) go to 50
    if(n48.lt.0) go to 50
    if(n44.lt.0) go to 50
    k     = n1 + n2 + n5 + n6 + n9 + n10 + n21 + n23 + n25 + n27 + &
          & n29 + n31 + n41 + n44 + n45 + n48 + n49 + n52
    if(k.gt.0) go to 54
    c9j  = 1.d0 
    return
50 c9j  = 0.d0
    return
!
!---    calcul de la somme sur xj
!
54 xn   = 2.d0*(xj21 - xj32) + eps1
    if(xn.lt.0.d0) xn  =  - (xn + eps2)
    jmin = xn
    xn   = 2.d0*(xj11 - xj33) + eps1
    if(xn.lt.0.d0) xn  =  - (xn + eps2)
    n     = xn
    if(n.gt.jmin) jmin = n
    xn   = 2.d0*(xj12 - xj23) + eps1
    if(xn.lt.0.d0) xn  =  - (xn + eps2)
    n     = xn
    if(n.gt.jmin) jmin = n
    xn   = 2.d0*(xj21 + xj32) + eps1
    jmax = xn
    xn   = 2.d0*(xj11 + xj33) + eps1
    n     = xn
    if(n.lt.jmax) jmax = n
    xn   = 2.d0*(xj12 + xj23) + eps1
    n     = xn
    if(n.lt.jmax) jmax = n
    xjmin  = dfloat(jmin)/2.d0
    xjmax  = dfloat(jmax)/2.d0
    s        = 0.d0
    xj   = xjmin
    xn   =  - xj32 + xj21 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n3   = xn
    xn   =  - xj11 + xj33 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n4   = xn
    xn   = xj32 - xj21 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n7   = xn
    xn   = xj11 - xj33 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n8   = xn
    xn   = xj32 + xj21 - xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n11  = xn
    xn   = xj11 + xj33 - xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n12  = xn
    xn   =  - xj11 - xj32 + xj31 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n13  = xn
    xn   =  - xj21 - xj33 + xj31 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n14  = xn
    xn   = xj32 + xj21 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n18  = xn
    n18  = n18 + 1
    xn   = xj11 + xj33 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n19  = xn
    n19  = n19 + 1
    xn   =  - xj21 + xj + xj32 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n22  = xn
    xn   =  - xj12 + xj23 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n24  = xn
    xn   = xj21 - xj + xj32 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n26  = xn
    xn   = xj12 - xj + xj23 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n28  = xn
    xn   = xj21 + xj - xj32 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n30  = xn
    xn   = xj12 + xj - xj23 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n32  = xn
    xn   =  - xj22 - xj + xj32 + xj23 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n34  = xn
    xn   = xj12 + xj21 + xj22 + xj + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n35  = xn
    n35  = n35 + 1
    xn   = xj21 + xj + xj32 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n37  = xn
    n37  = n37 + 1
    xn   = xj12 + xj + xj23 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n39  = xn
    n39  = n39 + 1
    xn   =  - xj + xj11 + xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n42  = xn
    xn   =  - xj + xj23 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n43  = xn
    xn   = xj - xj11 + xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n46  = xn
    xn   = xj - xj23 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n47  = xn
    xn   = xj + xj11 - xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n50  = xn
    xn   = xj + xj23 - xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n51  = xn
    xn   =  - xj13 - xj + xj33 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n53  = xn
    xn   = xj13 + xj + xj23 + xj11 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n55  = xn
    n55  = n55 + 1
    xn   = xj + xj11 + xj33 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n57  = xn
    n57  = n57 + 1
    xn   = xj + xj23 + xj12 + eps1
    if(xn.lt.0.d0) xn  = xn + eps2
    n58  = xn
    n58  = n58 + 1
    go to 10
    52 if(xj.gt.xjmax) go to 120
    n3   = n3 + 1
    n4   = n4 + 1
    n7   = n7 + 1
    n8   = n8 + 1
    n11  = n11 - 1
    n12  = n12 - 1
    n13  = n13 + 1
    n14  = n14 + 1
    n18  = n18 + 1
    n19  = n19 + 1
    n22  = n22 + 1
    n24  = n24 + 1
    n26  = n26 - 1
    n28  = n28 - 1
    n30  = n30 + 1
    n32  = n32 + 1
    n34  = n34 - 1
    n35  = n35 + 1
    n37  = n37 + 1
    n39  = n39 + 1
    n42  = n42 - 1
    n43  = n43 - 1
    n46  = n46 + 1
    n47  = n47 + 1
    n50  = n50 + 1
    n51  = n51 + 1
    n53  = n53 - 1
    n55  = n55 + 1
    n57  = n57 + 1
    n58  = n58 + 1
!
!--- calcul des sommes alternees s1(k1),s2(k2),s3(k3)
!
10 k1    = 0.d0
    l1   =  - n13
    if(l1.gt.k1)  k1    = l1
    l1   =  - n14
    if(l1.gt.k1)  k1    = l1
    l1   = n9
    if(n10.lt.l1) l1    = n10
    if(n11.lt.l1) l1    = n11
    if(n12.lt.l1) l1    = n12
    if(n15.lt.l1) l1    = n15
    f1   = 1.d0
    s1   = 1.d0
    i1   = k1 + 1
62 if(i1.gt.l1) go to 64
    i1m1     = i1 - 1
    nn1  = (n9 - i1m1)*(n10 - i1m1)*(n11 - i1m1)*(n12 - i1m1)
    nd1  = i1*(n13 + i1)*(n14 + i1)*(n15 - i1m1)
    f1   =  - f1*dfloat(nn1)/dfloat(nd1)
    s1   = s1 + f1
    i1   = i1 + 1
    go to 62
64 k2    = 0
    l2   =  - n33
    if(l2.gt.k2) k2  = l2
    l2   =  - n34
    if(l2.gt.k2) k2  = l2
    l2   = n29
    if(n30.lt.l2) l2    = n30
    if(n31.lt.l2) l2    = n31
    if(n32.lt.l2) l2    = n32
    if(n35.lt.l2) l2    = n35
    f2   = 1.d0
    s2   = 1.d0
    i2   = k2 + 1
70 if(i2.gt.l2) go to 80
    i2m2     = i2 - 1
    nn2  = (n29 - i2m2)*(n30 - i2m2)*(n31 - i2m2)*(n32 - i2m2)
    nd2  = i2*(n33 + i2)*(n34 + i2)*(n35 - i2m2)
    f2   =  - f2*dfloat(nn2)/dfloat(nd2)
    s2   = s2 + f2
    i2   = i2 + 1
    go to 70
80 k3    = 0
    l3   =  - n53
    if(l3.gt.k3) k3  = l3
    l3   =  - n54
    if(l3.gt.k3) k3  = l3
    l3   = n49
    if(n50.lt.l3) l3    = n50
    if(n51.lt.l3) l3    = n51
    if(n52.lt.l3) l3    = n52
    if(n55.lt.l3) l3    = n55
    f3   = 1.d0
    s3   = 1.d0
    i3   = k3 + 1
84 if(i3.gt.l3) go to 90
    i3m3     = i3 - 1
    nn3  = (n49 - i3m3)*(n50 - i3m3)*(n51 - i3m3)*(n52 - i3m3)
    nd3  = i3*(n53 + i3)*(n54 + i3)*(n55 - i3m3)
    f3   =  - f3*dfloat(nn3)/dfloat(nd3)
    s3   = s3 + f3
    i3   = i3 + 1
    go to 84
!
!--- calcul de la racine d'un terme de la somme sur j
!
90 s2n   = flog(n3+1)   + flog(n4+1)  + flog(n7+1)  + flog(n8+1)    + &
             & flog(n11+1) + flog(n12+1) + flog(n22+1) + flog(n24+1) + &
             & flog(n26+1) + flog(n28+1) + flog(n30+1) + flog(n32+1) + &
             & flog(n42+1) + flog(n43+1) + flog(n46+1) + flog(n47+1) + &
             & flog(n50+1) + flog(n51+1)
    s2n  = .5d0*s2n
    s2d  = flog(n18+1) + flog(n19+1) + flog(n37+1) + flog(n39+1) + &
             & flog(n57+1) + flog(n58+1)
    s2d  = .5d0*s2d
    km1  = k1 - 1
    kp1  = k1 + 1
    km2  = k2 - 1
    kp2  = k2 + 1
    km3  = k3 - 1
    kp3  = k3 + 1
    s2n  = s2n + flog(n15 - km1) + flog(n35 - km2) + flog(n55 - km3)
    s2d  = s2d + flog(kp1) + flog(kp2) + flog(kp3) + flog(n9 - km1) + &
             & flog(n10 - km1) + flog(n11 - km1) + flog(n12 - km1) + &
             & flog(n13 + kp1) + flog(n14 + kp1) + flog(n29 - km2) + &
             & flog(n30 - km2) + flog(n31 - km2) + flog(n32 - km2) + &
             & flog(n33 + kp2) + flog(n34 + kp2) + flog(n49 - km3) + &
             & flog(n50 - km3) + flog(n51 - km3) + flog(n52 - km3) + &
             & flog(n53 + kp3) + flog(n54 + kp3)
    f        = s2d - s2n
    if (f.gt.80.d0) go to 100
    f        = s2n/s2d
    if((f.lt.1.01d0).and.(f.gt.0.98d0)) go to 100
    f        = s1*s2*s3*dexp(s2n - s2d)
    go to 110
100 f    = s1*s2*s3
!    if(f) 102,112,104
    if(f.lt.0.0d0) then
        go to 102
    else if (f.eq.0.0d0) then
        go to 112
    else
        go to 104
    end if

102 f    = dlog( - f)
    f        =  - dexp(f + s2n - s2d)
    go to 110
104 f    = dlog(f)
    f        = dexp(f + s2n - s2d)
110 f    = f*(2.d0*xj + 1.d0)
!
!--- calcul de la phase d'un terme de la somme sur j
!
    l1   = k1 + k2 + k3
    k1   = l1/2
    k1   = 2*k1
    if(l1.ne.k1) f  =  - f
    s        = s + f
112 xj   = xj + 1.d0
    go to 52
!
!--- calcul du c9j sans phase   tro08270
!
120 c2n  = flog(n1+1) + flog(n2+1) + flog(n5+1) + flog(n6+1) + &
             & flog(n9+1) + flog(n10+1) + flog(n21+1) + flog(n23+1) + &
             & flog(n25+1) + flog(n27+1) + flog(n29+1) + flog(n31+1) + &
             & flog(n41+1) + flog(n44+1) + flog(n45+1) + flog(n48+1) + &
             & flog(n49+1) + flog(n52+1)
    c2n  = .5d0*c2n
    c2d  = flog(n16+1) + flog(n17+1) + flog(n36+1) + flog(n38+1) + &
             & flog(n56+1) + flog(n59+1)
    c2d  = .5d0*c2d
    f        = c2d - c2n
    if(f.gt.80.d0) go to 122
    f        = c2n/c2d
    if((f.lt.1.01d0).and.(f.gt.0.98d0)) go to 122
    c9j  = s*dexp(c2n - c2d)
    go to 130
!122 if(s) 124,50,126

122 if(s.lt.0.0d0) then
        go to 124
    else if(s.eq.0.0d0) then
        go to 50
    else
        go to 126
    end if

124 s    = dlog( - s)
    c9j  =  - dexp(s + c2n - c2d)
    go to 130
126 s    = dlog(s)
    c9j  = dexp(s + c2n - c2d)
!
!       calcul de la phase
!
130 k    = n9 + n16 + n36 + n56 - 1
    l        = k/2
    l        = 2*l
    if (l.ne.k) c9j  =  - c9j
    return
end subroutine neufj


!=============================================================================!
!
SUBROUTINE gfv
!
!=============================================================================!
!
!       Calculates sign, dsqrt, factorials, etc. of integers and half int.
!
!       iv(n)  =  (-1)**n
!       sq(n)  =  dsqrt(n)
!       sqi(n) =  1/dsqrt(n)
!       sqh(n) =  dsqrt(n+1/2)
!       shi(n) =  1/dsqrt(n+1/2)
!       fak(n) =  n!
!       fad(n) =  (2*n+1)!!
!       fi(n)  =  1/n!
!       wf(n)  =  dsqrt(n!)
!       wfi(n) =  1/dsqrt(n!)
!       gm2(n) =  gamma(n+1/2)
!       gmi(n) =  1/gamma(n+1/2)
!       wg(n)  =  dsqrt(gamma(n+1/2))
!       wgi(n) =  1/dsqrt(gamma(n+1/2))
!
!-----------------------------------------------------------------------
    
    integer :: i


!    third  = one/3.d0
!    pi   = 4*atan(one)

    iv(0)  = +1
    sq(0)  =  zero
    sqi(0) =  1.d30
    sqh(0) =  dsqrt(half)
    shi(0) =  1/sqh(0)
    fak(0) =  1
    fad(0) =  1
    fi(0)  =  1
    wf(0)  =  1
    wfi(0) =  1
!   gm2(0) = Gamma(1/2) = dsqrt(pi)
    gm2(0) =  dsqrt(pi)
    gmi(0) =  1/gm2(0)
    wg(0)  =  dsqrt(gm2(0))
    wgi(0) =  1/wg(0)
    do i = 1,igfv
        iv(i)  = -iv(i-1)
        sq(i)  = dsqrt(dfloat(i))
        sqi(i) = 1./sq(i)
        sqh(i) = dsqrt(i+half)
        shi(i) = 1./sqh(i)
        fak(i) = i*fak(i-1)
        fad(i) = (2*i+1.)*fad(i-1)
        fi(i)  = 1./fak(i)
        wf(i)  = sq(i)*wf(i-1)
        wfi(i) = 1./wf(i)
        gm2(i) = (i-half)*gm2(i-1)
        gmi(i) = 1./gm2(i)
        wg(i)  = sqh(i-1)*wg(i-1)
        wgi(i) = 1./wg(i)
    enddo
!
!       write(6,*) ' ****** END GFV *************************************' 
    return
!-end-GFV
end subroutine gfv

!=============================================================================!
!
FUNCTION racah0(j1,j2,j3,l1,l2,l3)
!
!=============================================================================!
!
!       Calculates 6j-symbol (notation of Edmonds) for integer ang.momenta
!
!-----------------------------------------------------------------------
    integer, intent(in) :: j1, j2, j3, l1, l2, l3
    

    double precision :: racah0

    integer :: i1, i2, i3, i4, i5, i6, i7, n, n1, n2
!
!    dsq(i,k,l) = wf(i+k-l)*wfi(i+k+l+1)*wf(i-k+l)*wf(k+l-i)
!
    racah0 = 0.0d0
    i1   = j1 + j2 + j3
    i2   = j1 + l2 + l3
    i3   = l1 + j2 + l3
    i4   = l1 + l2 + j3
    i5   = j1 + j2 + l1 + l2
    i6   = j2 + j3 + l2 + l3
    i7   = j3 + j1 + l3 + l1
    n1   = max0(i1,i2,i3,i4)
    n2   = min0(i5,i6,i7)

    if (n1.gt.n2) return
    do n = n1,n2
        racah0 = racah0  +  iv(n)*fak(n+1)*fi(n-i1)*fi(n-i2)*fi(n-i3)* &
                 & fi(n-i4)*fi(i5-n)*fi(i6-n)*fi(i7-n)
        
    end do
    racah0 = dsq(j1,j2,j3)*dsq(j1,l2,l3)*dsq(l1,j2,l3)*dsq(l1,l2,j3)*racah0
!
    return
! - end - RACAH0
end function racah0

!=============================================================================!
!
FUNCTION racslj(k,l1,l2,j2,j1)
!
!=============================================================================!
!
!       Calculates the Racah - coefficient      (   k  l1  l2 )
!                                               ( 1/2  j2  j1 )
!
!       for integer values               k = K, l1 = L1,    l2 = L2
!       and half integer values         j1 = J1 - 1/2, j2 = J2 - 1/2
!
!       Method of Edmonds
!

    double precision :: racslj
    integer, intent(in) :: k, l1, l2, j2, j1
    
    integer :: l12m, l12p, j12m, j12p

!
!    w(i,j,k,l,m,n) = sq(i)*sq(j)*sqi(k)*sqi(l)*sqi(m)*sqi(n) 
!
    racslj = 0.d0
!
    l12m     = l1  -  l2
    l12p     = l1  +  l2
    j12m     = j1  -  j2
    j12p     = j1  +  j2  -  1

!--- check of triangular rule
    if ( (iabs(l12m).gt.k .or. k.gt.l12p) .or. &
      &   (iabs(j12m).gt.k .or. k.gt.j12p) ) return
!
    if (j1.eq.l1 + 1) then
        if (j2.eq.l2 + 1) then
            racslj = -wsix(j12p+k+1, j12p-k,2*j1-1, j1, 2*j2-1, j2)
        elseif (j2.eq.l2) then
            racslj =  wsix(k-l12m, k+j12m, 2*j1-1, j1, l2, 2*l2+1)
        endif
    elseif (j1.eq.l1) then
        if (j2.eq.l2 + 1) then
            racslj =  wsix(k+l12m, k-j12m, 2*j2-1, j2, l1, 2*l1+1)
        elseif (j2.eq.l2) then
            racslj =  wsix(l12p+k+1, l12p-k, l1, 2*l1+1, l2, 2*l2+1)
        endif
    endif   
    racslj = iv(l12p+k)*racslj/2
!
    return
! - end - RACSLJ
end function racslj
!=============================================================================!
!
FUNCTION wiglll(l1,l2,l3)
!
!=============================================================================!
!
!       Calculates the Wigner - coefficient     ( l1    l2  l3 )
!                                              (   0    0   0 )
!       for integer values of l1,l2,l3
!
!       Method of Edmonds
!

    double precision :: wiglll
    integer,intent(in) :: l1, l2, l3

    integer :: l, l12p, l12m, lh
!
    wiglll = 0.d0
!
    l     = l1  +  l2   +   l3
    if(mod(l,2).ne.0) return
    l12p = l1  +  l2
    l12m = l1  -  l2
    lh   = l/2
!
!--- check of triangular rule
    if (iabs(l12m).gt.l3 .or. l3.gt.l12p) return
!
    wiglll = iv(lh)*wf(l12p-l3)*wf(l12m+l3)*wfi(l+1)*wf(l3-l12m)* &
      &     fak(lh)*fi(lh-l1)*fi(lh-l2)*fi(lh-l3)
!
    return
! - end - WIGLLL
end function wiglll

function dsq(i, k, l)
    double precision :: dsq
    integer, intent(in) :: i, k, l
    
    dsq = wf(i+k-l)*wfi(i+k+l+1)*wf(i-k+l)*wf(k+l-i)
end function dsq

function wsix(i, j, k,l, m, n)
    double precision :: wsix
    integer, intent(in) :: i, j, k, l, m, n

     wsix = sq(i)*sq(j)*sqi(k)*sqi(l)*sqi(m)*sqi(n) 
end function

!=============================================================================!
!
SUBROUTINE TROISJ(XJ1,XJ2,XJ3,XM1,XM2,XM3,C3J)
!
!=============================================================================!
!      IMPLICIT REAL*4 (A-H,O-Z)
!--- calculates the 3-j symbol: | xj1 xj2 xj3 |
!---                            | xm1 xm2 xm3 |
!      IMPLICIT REAL*8 (A-H,O-Z)
    DOUBLE PRECISION, DIMENSION(301) :: FLOG
    DOUBLE PRECISION, INTENT(IN) :: XJ1, XJ2, XJ3, XM1, XM2, XM3
    DOUBLE PRECISION, INTENT(OUT) :: C3J
!
!     MISE EN DATA DES LOG(FACTORIELLE) ET DE EPS
!
    DATA FLOG(2:31)/0.D0,.69314718D0,1.7917595D0,3.1780538D0,4.7874917D0, &
    &        6.5792511D0,8.5251613D0,10.604603D0,12.801827D0,15.104413D0, &
    &        17.502307D0,19.987214D0,22.552163D0,25.191221D0,27.899271D0, &
    &        30.671860D0,33.505072D0,36.395445D0,39.339884D0,42.335616D0, &
    &        45.380139D0,48.471180D0,51.606674D0,54.784729D0,58.003604D0, &
    &        61.261702D0,64.557537D0,67.889743D0,71.257038D0,74.658235D0/

    DATA FLOG(32:61)/78.092223D0,81.557959D0,85.054466D0,88.580827D0,92.136175D0, &
    &    95.719694D0,99.330612D0,102.96820D0,106.63176D0,110.32064D0,114.03421D0, &
    &    117.77188D0,121.53308D0,125.31727D0,129.12393D0,132.95257D0,136.80272D0, &
    &    140.67392D0,144.56574D0,148.47776D0,152.40959D0,156.36083D0,160.33112D0, &
    &    164.32011D0,168.32744D0,172.35279D0,176.39584D0,180.45629D0,184.53383D0, &
    &    188.62817D0/
    DATA FLOG(62:91)/192.73904D0,196.86618D0,201.00931D0,205.16820D0,209.34258D0, &
    &    213.53224D0,217.73693D0,221.95644D0,226.19054D0,230.43904D0,234.70172D0, &
    &    238.97839D0,243.26885D0,247.57291D0,251.89040D0,256.22113D0,260.56494D0, &
    &    264.92164D0,269.29110D0,273.67312D0,278.06757D0,282.47429D0,286.89313D0, &
    &    291.32394D0,295.76659D0,300.22094D0,304.68685D0,309.16419D0,313.65283D0, &
    &    318.15264D0/

    DATA FLOG(92:121)/322.66349D0,327.18529D0,331.71788D0,336.26118D0,340.81505D0,&
    &     345.37940D0,349.95411D0,354.53908D0,359.13420D0,363.73937D0,368.35449D0,&
    &     372.97946D0,377.61419D0,382.25859D0,386.91255D0,391.57598D0,396.24881D0,&
    &     400.93094D0,405.62230D0,410.32277D0,415.03230D0,419.75080D0,424.47819D0,&
    &     429.21439D0,433.95932D0,438.71291D0,443.47508D0,448.24576D0,453.02489D0,&
    &     457.81238D0/

    DATA FLOG(122:151)/462.60817D0,467.41220D0,472.22438D0,477.04466D0,481.87298D0,&
    &      486.70926D0,491.55345D0,496.40547D0,501.26529D0,506.13282D0,511.00802D0,&
    &      515.89082D0,520.78117D0,525.67901D0,530.58428D0,535.49694D0,540.41692D0,&
    &      545.34417D0,550.27865D0,555.22029D0,560.16905D0,565.12488D0,570.08772D0,&
    &      575.05753D0,580.03427D0,585.01787D0,590.00830D0,595.00552D0,600.00946D0,&
    &      605.02010D0/

    DATA FLOG(152:181)/610.03738D0,615.06126D0,620.09170D0,625.12866D0,630.17208D0,&
    &      635.22193D0,640.27818D0,645.34077D0,650.40968D0,655.48486D0,660.56626D0,&
    &      665.65385D0,670.74760D0,675.84747D0,680.95341D0,686.06541D0,691.18340D0,&
    &      696.30735D0,701.43726D0,706.57306D0,711.71472D0,716.86221D0,722.01551D0,&
    &      727.17456D0,732.33934D0,737.50983D0,742.68598D0,747.86776D0,753.05516D0,&
    &      758.24811D0/

    DATA FLOG(182:211)/763.44661D0,768.65061D0,773.86010D0,779.07503D0,784.29539D0,&
    &      789.52114D0,794.75224D0,799.98869D0,805.23044D0,810.47747D0,815.72973D0,&
    &      820.98722D0,826.24991D0,831.51778D0,836.79078D0,842.06890D0,847.35209D0,&
    &      852.64036D0,857.93366D0,863.23199D0,868.53529D0,873.84356D0,879.15676D0,&
    &      884.47488D0,889.79789D0,895.12577D0,900.45848D0,905.79603D0,911.13836D0,&
    &      916.48547D0/

    DATA FLOG(212:241)/921.83732D0,927.19391D0,932.55521D0,937.92118D0,943.29181D0,&
    &      948.66710D0,954.04699D0,959.43148D0,964.82056D0,970.21419D0,975.61235D0,&
    &      981.01503D0,986.42220D0,991.83385D0,997.24995D0,1002.6705D0,1008.0954D0,&
    &      1013.5248D0,1018.9585D0,1024.3966D0,1029.8389D0,1035.2857D0,1040.7367D0,&
    &      1046.1920D0,1051.6516D0,1057.1155D0,1062.5836D0,1068.0558D0,1073.5323D0,&
    &      1079.0129D0/

    DATA FLOG(242:271)/1084.4977D0,1089.9866D0,1095.4797D0,1100.9768D0,1106.4781D0,&
    &      1111.9834D0,1117.4928D0,1123.0063D0,1128.5237D0,1134.0452D0,1139.5706D0,&
    &      1145.1001D0,1150.6335D0,1156.1708D0,1161.7120D0,1167.2573D0,1172.8063D0,&
    &      1178.3593D0,1183.9161D0,1189.4768D0,1195.0413D0,1200.6097D0,1206.1818D0,&
    &      1211.7577D0,1217.3375D0,1222.9209D0,1228.5082D0,1234.0992D0,1239.6939D0,&
    &      1245.2924D0/
    
    DATA FLOG(272:301)/1250.8944D0,1256.5003D0,1262.1097D0,1267.7228D0,1273.3396D0,&
    &      1278.9600D0,1284.5840D0,1290.2117D0,1295.8429D0,1301.4777D0,1307.1160D0,&
    &      1312.7580D0,1318.4034D0,1324.0524D0,1329.7048D0,1335.3609D0,1341.0203D0,&
    &      1346.6833D0,1352.3497D0,1358.0196D0,1363.6929D0,1369.3697D0,1375.0499D0,&
    &      1380.7334D0,1386.4204D0,1392.1107D0,1397.8045D0,1403.5016D0,1409.2020D0,&
    &      1414.9058D0/

    DOUBLE PRECISION :: EPS1, EPS2, XN, xiii, F, S, C2N, C2D
    INTEGER :: N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, K, L, I, IM1, NN
    INTEGER :: KM1, KP1, ND

    DATA EPS1,EPS2/.1D0,-.2D0/
!
!     CALCUL DES COMBINAISONS J,M
!
    XN  =  XJ2-XM2+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N1  =  XN
    xiii=  dabs(xj1)+dabs(xj2)+dabs(xj3)
    if(xiii.lt.200) go to 9999
    write (11,*) 'xj1 xj2 xj3=  ', xj1,xj2,xj3
    stop 
 9999 continue
    XN  =  XJ3+XM3+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N2  =  XN
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    XN  =  XJ3-XM3+EPS1
    N3  =  XN
    XN  =  XJ1+XM1+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N4  =  XN
    XN  =  XJ2+XJ3-XJ1+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N5  =  XN
    XN  =  XJ1+XJ3-XJ2+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N6  =  XN
    XN  =  XJ1+XJ2-XJ3+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N7  =  XN
    XN  =  XJ1-XM1+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N8  =  XN
    XN  =  XJ2+XM2+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N9  =  XN
    XN  =  XJ3-XJ1-XM2+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N10 =  XN
    XN  =  XJ3-XJ2+XM1+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N11 =  XN
    XN  =  XJ1+XJ2+XJ3+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    N12 =  XN
    N12 =  N12+1
!
!     TESTS SUR LES J ET M
!
    K   =  N4*N8
    IF(K.LT.0) GO TO 59
    K   =  N1*N9
    IF(K.LT.0) GO TO 59
    K   =  N2*N3
    IF(K.LT.0) GO TO 59
    IF(N5.LT.0) GO TO 59
    IF(N6.LT.0) GO TO 59
    IF(N7.LT.0) GO TO 59
    L   =  N1-N2+N3-N4+N8-N9
    IF(L.NE.0) GO TO 59
    K   =  N12-1
    IF(K.GT.0) GO TO 60
    C3J =  1.D0
    RETURN
 59 C3J =  0.D0
    RETURN
!
!     CALCUL DE LA SOMME ALTERNEE
!
 60 K   =  0
    L   =  -N10
    IF(L.GT.K) K    =  L
    L   =  -N11
    IF(L.GT.K) K    =  L
    L   =  N7
    IF(N8.GT.L) L   =  N8
    IF(N9.GT.L) L   =  N9
    F   =  1.D0
    S   =  1.D0
    I   =  K+1
   62 IF(I.GT.L) GO TO 80
    IM1 =  I-1
    NN  =  (N7-IM1)*(N8-IM1)*(N9-IM1)
    ND  =  I*(N10+I)*(N11+I)
    F   =  -F*DFLOAT(NN)/DFLOAT(ND)
    S   =  S+F
    I   =  I+1
    GO TO 62
!
!     CALCUL DE LA RACINE
!
 80 C2N =  FLOG(N1+1) + FLOG(N2+1) + FLOG(N3+1) + FLOG(N4+1) + FLOG(N5+1) + FLOG(N6+1) + &
    &      FLOG(N7+1) + FLOG(N8+1) + FLOG(N9+1)
    C2N =  .5D0*C2N
    KM1 =  K-1
    KP1 =  K+1
    C2D =  FLOG(KP1) + FLOG(N7-KM1) + FLOG(N8-KM1) + FLOG(N9-KM1) + FLOG(N10+KP1) + &
    &      FLOG(N11+KP1) + .5D0*FLOG(N12+1)
!
!     CALCUL DU C3J SANS PHASE
!
    F   =  C2D-C2N
    IF(F.GT.80.D0) GO TO 98
    F   =  C2N/C2D
    IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 98
    C3J =  S*DEXP(C2N-C2D)
    GO TO 106
!   98 IF(S) 100,59,102
 98 IF(S.lt.0.0d0) then
        go to 100
     else if(S.eq.0.0d0) then
        go to 59
     else
        go to 102
    end if

100 S   =  DLOG(-S)
    C3J =  -DEXP(S+C2N-C2D)
    GO TO 106
102 S   =  DLOG(S)
    C3J =  DEXP(S+C2N-C2D)
!
!     CALCUL DE LA PHASE
!
106 XN  =  XJ1-XJ2-XM3+EPS1
    IF(XN.LT.0.D0) XN   =  XN+EPS2
    L   =  XN
    L   =  L+K
    K   =  L/2
    K   =  2*K
    IF(L.NE.K) C3J  =  -C3J
    RETURN
END SUBROUTINE TROISJ


end module jsymbols