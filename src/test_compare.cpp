#include "classgroup.h"
#include "classgroupSquarer.h"
#include "timer.h"
#include <iostream>

int main(int argc, char* argv[]) {
  const int num_iters = 10000;
  timer time_calc;

  std::cout << "benchmark small 512 bit prime\n";
  // relic_int oldD = -77447;
  relic_int D;

  for (int i = 0; i < 3; i++) {
    switch (i) {
      case 0:
        D = "-117486081361718979337668786659568125756677310916264821465703349385729177030045693203384244824767784173284"
            "73277751253005897383252417809294541920707334384671";
        break;
      case 1:
        D = "-157523081333953688776525695336976710923610890627472321818455934829777889400307231472509445273569926623020"
            "9155591213712581215593075380505087533130123281860532174327789052159663080760147218557354963529169613391033"
            "17220859109143821373938156368140794887760611548625948124400804378932790799247109454639983606744623";
        break;
      case 2:
        D = "-283178161896406766256338548801118480561159510501921419471614013957318883445267258354383871905867382312137"
            "1050638064399792114486970359147566019220785055992688075893732299055316117066052956001296964804915301153886"
            "7036290506786163373332067953476686329787104129027881326169410875103254720158321671870574669799168483168258"
            "5497782293394919544553844305624763702105253463746988655503565742698689196724684751881864328929476448574800"
            "4297312897711426338328349889340902253662139880141610848900485391388016368326796006101572522791522472520876"
            "4732811757663881381962068567438737166388615608125137092107866274370573147791660110537263";
        break;
    }

    std::cout << "D = " << D << " " << D.bitSize() << "\n";

    relic_int one(1);
    relic_int two(2);
    Classgroup g = from_ab_discriminant(two, one, D);
    // std::cout << "Starting g = " << g << " {" << g.a_bits() << "," << g.b_bits() << "," << g.c_bits() << "}\n";

    while ((g.a_bits() < g.c_bits()) || (g.b_bits() < g.c_bits())) { g = pow(g, 2); }

    std::cout << "g = " << g << "\n";
    Classgroup g2 = pow(g, 2);

    std::cout << "------------------------------------------------Start comparison......................\n";

    ClassgroupSquarer Sqr(g2);

    time_calc.start();
    for (int i = 0; i < num_iters; i++) { g2 = square(g2); }
    time_calc.stop();
    std::cout << "Time Elapsed for square(g2): " << time_calc.duration_in_milliseconds(num_iters)
              << " ms per iteration\t" << time_calc.duration_in_milliseconds() << " ms total \n";

    time_calc.start();
    for (int i = 0; i < num_iters; i++) { Sqr.square(); }
    time_calc.stop();
    std::cout << "Time Elapsed for Sqr.square(): " << time_calc.duration_in_milliseconds(num_iters)
              << " ms per iteration\t" << time_calc.duration_in_milliseconds() << " ms total \n";

    Classgroup gSqr = Sqr.result();

    std::cout
        << "Final ------------------------------------------------------------------------------------------------\n";
    std::cout << "Ref g2 = " << g2 << "\n";
    std::cout << "------------------------------------------------------------------------------------------------\n";
    std::cout << "New g2 = " << gSqr << "\n";
    std::cout
        << "Done ------------------------------------------------------------------------------------------------\n";

    if (g2 == gSqr)
      std::cout << "\n------------MATCH-------------\n";
    else
      std::cout << "\n------------FAIL-------------\n";

    //--------------------------------------------------------------------------------------------------
    double sum_time = 0;
    for (int j = 0; j < 10; j++) {
      time_calc.start();
      for (int i = 0; i < num_iters; i++) Sqr.square();
      time_calc.stop();
      sum_time += time_calc.duration_in_milliseconds();
      std::cout << "Time Elapsed for Sqr.square(): " << time_calc.duration_in_milliseconds(num_iters)
                << " ms per iteration\t" << time_calc.duration_in_milliseconds() << " ms total \n";
    }

    gSqr = Sqr.result();

    double Sqr_time = sum_time;
    std::cout << "Sum time = " << sum_time / 1000.0 << "s \n";
    std::cout << "Final Sqr g2 =  " << gSqr << "\n";

    //--------------------------------------------------------------------------------------------------
    sum_time = 0;
    for (int j = 0; j < 10; j++) {
      time_calc.start();
      for (int i = 0; i < num_iters; i++) g2 = square(g2);
      time_calc.stop();
      sum_time += time_calc.duration_in_milliseconds();
      std::cout << "Time Elapsed for g.square(): " << time_calc.duration_in_milliseconds(num_iters)
                << " ms per iteration\t" << time_calc.duration_in_milliseconds() << " ms total \n";
    }
    std::cout << "Sum time = " << sum_time / 1000.0 << "s \n";
    std::cout << "Final g2 =  " << g2 << "\n";
    //--------------------------------------------------------------------------------------------------

    double ratio = int(1000.0 * Sqr_time / sum_time) / 10.0;

    if (g2 == gSqr)
      std::cout << "\n------------MATCH-------------\n";
    else {
      std::cout << "\n------------FAIL-------------\n";
      return 0;
    }

    std::cout << "********************* Speed factor = " << ratio << "%************************\n";

    //---------------------------------------------------------------------------------------------------
    /*
    if (i == 0) {
      relic_int a("25596790638743725958110290155156497222819653591454105219587905336875299016469");
      relic_int b("-23832021859816254426257925958364506181090204426684233226584511847595213095497");
      relic_int c("120294117883045153967207226168660553442119760625332200739014312594578680079180");
      Classgroup final_g2(a, b, c);

      if (g2 == final_g2)
        std::cout << " matches expected result\n";
      else
        std::cout << " ERROR *********** failed to match \n";
    }
    */
  }
  std::cout << "\n------------SUCCESS-------------\n";
  return 0;
}
