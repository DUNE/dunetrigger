
#include "dunetrigger/triggeralgs/include/triggeralgs/BSMWindow/models/treelite_compmodel_classifier_xgboost/treelitemodel.h"


static const int32_t num_class[] = {  1, };

int32_t get_num_target(void) {
  return N_TARGET;
}
void get_num_class(int32_t* out) {
  for (int i = 0; i < N_TARGET; ++i) {
    out[i] = num_class[i];
  }
}
int32_t get_num_feature(void) {
  return 25;
}
const char* get_threshold_type(void) {
  return "float32";
}
const char* get_leaf_output_type(void) {
  return "float32";
}

void predict(union Entry* data, int pred_margin, float* result) {
  unsigned int tmp;
  if ( (data[22].missing != -1) && (data[22].fvalue < (float)154247)) {
    if ( (data[13].missing != -1) && (data[13].fvalue < (float)342922)) {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)692)) {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)126631)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)151973)) {
            result[0] += -0.6001512;
          } else {
            result[0] += -0.21829484;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)212772)) {
            result[0] += 0.3115297;
          } else {
            result[0] += -0.26979515;
          }
        }
      } else {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)164153)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)364428)) {
            result[0] += -0.26981968;
          } else {
            result[0] += 0.45023656;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)33528)) {
            result[0] += 0.4233306;
          } else {
            result[0] += -0.11439011;
          }
        }
      }
    } else {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)127)) {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)761883)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)526236)) {
            result[0] += -0.42530242;
          } else {
            result[0] += 0.08427634;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)251792)) {
            result[0] += -0.11829113;
          } else {
            result[0] += 0.7631764;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)626)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)673688)) {
            result[0] += -0.48100233;
          } else {
            result[0] += 0.67057747;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)684)) {
            result[0] += 0.25483173;
          } else {
            result[0] += 0.68134415;
          }
        }
      }
    }
  } else {
    if ( (data[17].missing != -1) && (data[17].fvalue < (float)195846)) {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)410077)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)469633)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)91695)) {
            result[0] += -0.01957902;
          } else {
            result[0] += -0.54443395;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)638969)) {
            result[0] += 0.34133598;
          } else {
            result[0] += 0.72392434;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)411)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)986)) {
            result[0] += 0.21459809;
          } else {
            result[0] += -0.7418182;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)102285)) {
            result[0] += 0.6953905;
          } else {
            result[0] += -0.012826927;
          }
        }
      }
    } else {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)64)) {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)164265)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)1347)) {
            result[0] += -0.5213414;
          } else {
            result[0] += 0.21392691;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)220466)) {
            result[0] += 0.5918153;
          } else {
            result[0] += -0.6117647;
          }
        }
      } else {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)189)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)84871)) {
            result[0] += 0.51438534;
          } else {
            result[0] += -0.10095153;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)528993)) {
            result[0] += 0.665643;
          } else {
            result[0] += 0.78454596;
          }
        }
      }
    }
  }
  if ( (data[22].missing != -1) && (data[22].fvalue < (float)174012)) {
    if ( (data[8].missing != -1) && (data[8].fvalue < (float)392049)) {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)151886)) {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)190877)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)150997)) {
            result[0] += -0.38533077;
          } else {
            result[0] += 0.109146915;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)98975)) {
            result[0] += -0.23383875;
          } else {
            result[0] += 0.36783868;
          }
        }
      } else {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)70939)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)39954)) {
            result[0] += 0.3117437;
          } else {
            result[0] += -0.05023075;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)402866)) {
            result[0] += -0.31622344;
          } else {
            result[0] += 0.39579278;
          }
        }
      }
    } else {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)3403)) {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)147035)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)642320)) {
            result[0] += -0.6712675;
          } else {
            result[0] += -0.28412262;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)174005)) {
            result[0] += 0.53515595;
          } else {
            result[0] += -0.679392;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)214)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)70939)) {
            result[0] += 0.25400823;
          } else {
            result[0] += -0.49989754;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)477)) {
            result[0] += 0.17781061;
          } else {
            result[0] += 0.5889107;
          }
        }
      }
    }
  } else {
    if ( (data[22].missing != -1) && (data[22].fvalue < (float)507262)) {
      if ( (data[18].missing != -1) && (data[18].fvalue < (float)444532)) {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)182956)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)81701)) {
            result[0] += 0.12957752;
          } else {
            result[0] += -0.3876753;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)505667)) {
            result[0] += 0.25341603;
          } else {
            result[0] += 0.5313121;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)207)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)20274)) {
            result[0] += -0.76185465;
          } else {
            result[0] += 0.17486914;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)472)) {
            result[0] += 0.22829525;
          } else {
            result[0] += 0.577362;
          }
        }
      }
    } else {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)315)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)36696)) {
          result[0] += -0.75451607;
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)60194)) {
            result[0] += 0.33314043;
          } else {
            result[0] += -0.6139087;
          }
        }
      } else {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)693702)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)63)) {
            result[0] += -0.43488455;
          } else {
            result[0] += 0.49673843;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)522)) {
            result[0] += 0.12809025;
          } else {
            result[0] += 0.5765118;
          }
        }
      }
    }
  }
  if ( (data[18].missing != -1) && (data[18].fvalue < (float)354754)) {
    if ( (data[22].missing != -1) && (data[22].fvalue < (float)311373)) {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)255869)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)329)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)386866)) {
            result[0] += -0.314727;
          } else {
            result[0] += 0.1664781;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)30220)) {
            result[0] += 0.14197366;
          } else {
            result[0] += -0.19596112;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)1178)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)858395)) {
            result[0] += -0.24273889;
          } else {
            result[0] += 0.61503345;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)875)) {
            result[0] += 0.144118;
          } else {
            result[0] += 0.42683944;
          }
        }
      }
    } else {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)129)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)17437)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)733)) {
            result[0] += -0.22678979;
          } else {
            result[0] += 0.37891042;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)930)) {
            result[0] += -0.31054345;
          } else {
            result[0] += -0.81882447;
          }
        }
      } else {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)25658)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)138084)) {
            result[0] += 0.44345376;
          } else {
            result[0] += 0.17254294;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)693702)) {
            result[0] += -0.06653426;
          } else {
            result[0] += 0.40185758;
          }
        }
      }
    }
  } else {
    if ( (data[23].missing != -1) && (data[23].fvalue < (float)223)) {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)72363)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)177475)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)423581)) {
            result[0] += -0.53268725;
          } else {
            result[0] += 0.14161482;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)769)) {
            result[0] += 0.5763537;
          } else {
            result[0] += -0.5603989;
          }
        }
      } else {
        result[0] += -0.6755764;
      }
    } else {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)644)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)32777)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)62)) {
            result[0] += -0.65045416;
          } else {
            result[0] += 0.23679177;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)588030)) {
            result[0] += -0.7526797;
          } else {
            result[0] += 0.27496296;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)282)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)121103)) {
            result[0] += 0.21859427;
          } else {
            result[0] += -0.73213714;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)199)) {
            result[0] += 0.2014926;
          } else {
            result[0] += 0.49200168;
          }
        }
      }
    }
  }
  if ( (data[20].missing != -1) && (data[20].fvalue < (float)692)) {
    if ( (data[12].missing != -1) && (data[12].fvalue < (float)654)) {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)351260)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)827244)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)281)) {
            result[0] += -0.38526753;
          } else {
            result[0] += -0.14109325;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)46068)) {
            result[0] += -0.48728737;
          } else {
            result[0] += 0.5889546;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)192)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)130171)) {
            result[0] += -0.35200235;
          } else {
            result[0] += 0.5447056;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)1153)) {
            result[0] += -0.2758721;
          } else {
            result[0] += 0.6124675;
          }
        }
      }
    } else {
      if ( (data[7].missing != -1) && (data[7].fvalue < (float)179804)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)42335)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)563)) {
            result[0] += -0.17532355;
          } else {
            result[0] += 0.22365916;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)63425)) {
            result[0] += -0.13353032;
          } else {
            result[0] += -0.3862985;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)522)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)437876)) {
            result[0] += -0.005003996;
          } else {
            result[0] += 0.54653585;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)139489)) {
            result[0] += 0.5046672;
          } else {
            result[0] += -0.22233884;
          }
        }
      }
    }
  } else {
    if ( (data[17].missing != -1) && (data[17].fvalue < (float)367906)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)52047)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)347883)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)49256)) {
            result[0] += 0.16812943;
          } else {
            result[0] += -0.089749984;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)781)) {
            result[0] += -0.29170334;
          } else {
            result[0] += 0.43352276;
          }
        }
      } else {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)479164)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)77094)) {
            result[0] += -0.09901879;
          } else {
            result[0] += -0.413323;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)831)) {
            result[0] += -0.3958801;
          } else {
            result[0] += 0.5082059;
          }
        }
      }
    } else {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)296)) {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)28443)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)2339)) {
            result[0] += -0.19605392;
          } else {
            result[0] += 0.5498337;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)30840)) {
            result[0] += -0.62300587;
          } else {
            result[0] += 0.14439623;
          }
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)644)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)4114)) {
            result[0] += 0.18472749;
          } else {
            result[0] += -1.1761864;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)65)) {
            result[0] += 0.11153245;
          } else {
            result[0] += 0.4352761;
          }
        }
      }
    }
  }
  if ( (data[12].missing != -1) && (data[12].fvalue < (float)444613)) {
    if ( (data[23].missing != -1) && (data[23].fvalue < (float)517490)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)50251)) {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)56063)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)494)) {
            result[0] += -0.025349293;
          } else {
            result[0] += 0.20574498;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)356962)) {
            result[0] += -0.23047948;
          } else {
            result[0] += 0.33101046;
          }
        }
      } else {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)378416)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)579422)) {
            result[0] += -0.2760261;
          } else {
            result[0] += 0.47627336;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)930)) {
            result[0] += -0.5064645;
          } else {
            result[0] += 0.34370893;
          }
        }
      }
    } else {
      if ( (data[5].missing != -1) && (data[5].fvalue < (float)226218)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)1404)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)280876)) {
            result[0] += -0.6770659;
          } else {
            result[0] += 0.033331733;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)377421)) {
            result[0] += 0.41421542;
          } else {
            result[0] += -0.9045333;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)259163)) {
          result[0] += -0.99774736;
        } else {
          result[0] += 0.24010439;
        }
      }
    }
  } else {
    if ( (data[15].missing != -1) && (data[15].fvalue < (float)1519)) {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)81462)) {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)5196)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)190262)) {
            result[0] += 0.44495377;
          } else {
            result[0] += -0.6727318;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)409248)) {
            result[0] += -0.6286432;
          } else {
            result[0] += 0.20714748;
          }
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)6547)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)198818)) {
            result[0] += -0.55657697;
          } else {
            result[0] += -1.081701;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)21192)) {
            result[0] += -0.73289156;
          } else {
            result[0] += 0.22213577;
          }
        }
      }
    } else {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)309)) {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)43637)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)444)) {
            result[0] += -0.5740082;
          } else {
            result[0] += 0.25079867;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)195673)) {
            result[0] += -0.7420783;
          } else {
            result[0] += -0.08290399;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)128)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)148804)) {
            result[0] += 0.25496128;
          } else {
            result[0] += -0.46147177;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)193)) {
            result[0] += -0.27961332;
          } else {
            result[0] += 0.450773;
          }
        }
      }
    }
  }
  if ( (data[17].missing != -1) && (data[17].fvalue < (float)528993)) {
    if ( (data[7].missing != -1) && (data[7].fvalue < (float)520531)) {
      if ( (data[9].missing != -1) && (data[9].fvalue < (float)21609)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)39954)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)463)) {
            result[0] += -0.0033141428;
          } else {
            result[0] += 0.2320025;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)342922)) {
            result[0] += -0.22888508;
          } else {
            result[0] += 0.30889982;
          }
        }
      } else {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)601313)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)79383)) {
            result[0] += -0.12648721;
          } else {
            result[0] += -0.35719842;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)3623)) {
            result[0] += -0.2484703;
          } else {
            result[0] += 0.42063054;
          }
        }
      }
    } else {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)488)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)99492)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)144953)) {
            result[0] += -0.6998182;
          } else {
            result[0] += 0.065502144;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)377829)) {
            result[0] += -0.40676036;
          } else {
            result[0] += 0.5449916;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)1009)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)9768)) {
            result[0] += -0.45795974;
          } else {
            result[0] += 0.2927311;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)293)) {
            result[0] += -0.6217403;
          } else {
            result[0] += 0.44006452;
          }
        }
      }
    }
  } else {
    if ( (data[23].missing != -1) && (data[23].fvalue < (float)392)) {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)376389)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)90243)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)128)) {
            result[0] += 0.2572964;
          } else {
            result[0] += -0.42592245;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)264099)) {
            result[0] += -0.8366441;
          } else {
            result[0] += -0.35835576;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)85736)) {
          result[0] += -0.30423662;
        } else {
          result[0] += 0.5367986;
        }
      }
    } else {
      if ( (data[18].missing != -1) && (data[18].fvalue < (float)1986)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)61)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)392)) {
            result[0] += 0.39548758;
          } else {
            result[0] += -0.46355873;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)1261366)) {
            result[0] += -0.902232;
          } else {
            result[0] += 0.22114158;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)499)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)197078)) {
            result[0] += 0.25930092;
          } else {
            result[0] += -0.6206755;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)61)) {
            result[0] += 0.079628475;
          } else {
            result[0] += 0.4086791;
          }
        }
      }
    }
  }
  if ( (data[17].missing != -1) && (data[17].fvalue < (float)128)) {
    if ( (data[16].missing != -1) && (data[16].fvalue < (float)446)) {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)296814)) {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)387390)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)22565)) {
            result[0] += -0.16153315;
          } else {
            result[0] += -0.35471064;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)128)) {
            result[0] += 0.58861154;
          } else {
            result[0] += -0.5693125;
          }
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)93204)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)243391)) {
            result[0] += 0.32422814;
          } else {
            result[0] += -0.5946452;
          }
        } else {
          result[0] += -0.79316926;
        }
      }
    } else {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)692)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)25910)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)57370)) {
            result[0] += -0.08110127;
          } else {
            result[0] += 0.48312184;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)84947)) {
            result[0] += -0.45027953;
          } else {
            result[0] += -0.1283757;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)10226)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)222305)) {
            result[0] += 0.26712888;
          } else {
            result[0] += -0.25147876;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)61992)) {
            result[0] += 0.18696904;
          } else {
            result[0] += -0.30873865;
          }
        }
      }
    }
  } else {
    if ( (data[19].missing != -1) && (data[19].fvalue < (float)42465)) {
      if ( (data[14].missing != -1) && (data[14].fvalue < (float)67938)) {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)68756)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)132)) {
            result[0] += 0.13983186;
          } else {
            result[0] += 0.3924285;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)72912)) {
            result[0] += 0.19182365;
          } else {
            result[0] += -0.13777106;
          }
        }
      } else {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)561)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)40919)) {
            result[0] += -0.13367191;
          } else {
            result[0] += -0.53613514;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)132247)) {
            result[0] += 0.03293696;
          } else {
            result[0] += 0.4186276;
          }
        }
      }
    } else {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)199159)) {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)556445)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)27287)) {
            result[0] += 0.17134233;
          } else {
            result[0] += -0.25769007;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)1833)) {
            result[0] += -0.33110413;
          } else {
            result[0] += 0.35319278;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)336)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)7379)) {
            result[0] += -0.38959345;
          } else {
            result[0] += 0.003710263;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)61)) {
            result[0] += -0.13697287;
          } else {
            result[0] += 0.35854936;
          }
        }
      }
    }
  }
  if ( (data[22].missing != -1) && (data[22].fvalue < (float)693702)) {
    if ( (data[4].missing != -1) && (data[4].fvalue < (float)101789)) {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)386)) {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)94033)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)183435)) {
            result[0] += -0.019972073;
          } else {
            result[0] += -0.35381603;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)268225)) {
            result[0] += -0.33238953;
          } else {
            result[0] += 0.15330608;
          }
        }
      } else {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)30220)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)498)) {
            result[0] += 0.026251009;
          } else {
            result[0] += 0.2940985;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)211)) {
            result[0] += -0.17636472;
          } else {
            result[0] += 0.059203353;
          }
        }
      }
    } else {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)334070)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)389945)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)64663)) {
            result[0] += -0.21916442;
          } else {
            result[0] += -0.46250802;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)2080)) {
            result[0] += -0.45720148;
          } else {
            result[0] += 0.39868456;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)295)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)451043)) {
            result[0] += -0.558318;
          } else {
            result[0] += -0.08928733;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)7961)) {
            result[0] += -0.252701;
          } else {
            result[0] += 0.33638003;
          }
        }
      }
    }
  } else {
    if ( (data[23].missing != -1) && (data[23].fvalue < (float)4716)) {
      if ( (data[24].missing != -1) && (data[24].fvalue < (float)122)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)1442393)) {
          result[0] += -1.0912927;
        } else {
          result[0] += -0.092340685;
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)6078)) {
          result[0] += 0.37965044;
        } else {
          result[0] += -0.31180963;
        }
      }
    } else {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)41495)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)343219)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)166765)) {
            result[0] += 0.4231855;
          } else {
            result[0] += 0.18011884;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)367906)) {
            result[0] += -0.4762709;
          } else {
            result[0] += 0.3997564;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)371)) {
          result[0] += -0.52419585;
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)309)) {
            result[0] += 0.08405673;
          } else {
            result[0] += 0.37939256;
          }
        }
      }
    }
  }
  if ( (data[23].missing != -1) && (data[23].fvalue < (float)469633)) {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)691562)) {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)125)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)180934)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)163586)) {
            result[0] += -0.19718541;
          } else {
            result[0] += 0.050918277;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)99133)) {
            result[0] += 0.116722584;
          } else {
            result[0] += -0.40699005;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)36668)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)656)) {
            result[0] += 0.038010634;
          } else {
            result[0] += 0.2938231;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)559356)) {
            result[0] += -0.10992127;
          } else {
            result[0] += 0.2929604;
          }
        }
      }
    } else {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)376)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)94644)) {
          result[0] += -0.538247;
        } else {
          result[0] += 0.44020602;
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)48107)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)347087)) {
            result[0] += 0.6284283;
          } else {
            result[0] += 0.39551446;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)11562)) {
            result[0] += 0.38915148;
          } else {
            result[0] += -0.1917709;
          }
        }
      }
    }
  } else {
    if ( (data[20].missing != -1) && (data[20].fvalue < (float)131)) {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)909421)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)1096)) {
          result[0] += -0.72058386;
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)7078)) {
            result[0] += 0.078166954;
          } else {
            result[0] += -0.61219275;
          }
        }
      } else {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)49186)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)1296218)) {
            result[0] += 0.48495778;
          } else {
            result[0] += -0.2814413;
          }
        } else {
          result[0] += -0.4936052;
        }
      }
    } else {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)909421)) {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)260466)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)150224)) {
            result[0] += 0.28420573;
          } else {
            result[0] += -0.28429112;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)527264)) {
            result[0] += -0.83215123;
          } else {
            result[0] += 0.2583091;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)10969)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)303)) {
            result[0] += -0.6291624;
          } else {
            result[0] += 0.382035;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)211959)) {
            result[0] += 0.4412888;
          } else {
            result[0] += 0.26203096;
          }
        }
      }
    }
  }
  if ( (data[4].missing != -1) && (data[4].fvalue < (float)84822)) {
    if ( (data[9].missing != -1) && (data[9].fvalue < (float)84106)) {
      if ( (data[14].missing != -1) && (data[14].fvalue < (float)100013)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)132)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)84408)) {
            result[0] += 0.055587996;
          } else {
            result[0] += -0.279211;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)1183)) {
            result[0] += 0.045386743;
          } else {
            result[0] += 0.22153674;
          }
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)526236)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)38654)) {
            result[0] += 0.08531276;
          } else {
            result[0] += -0.22231548;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)1217)) {
            result[0] += 0.0020240725;
          } else {
            result[0] += 0.3939981;
          }
        }
      }
    } else {
      if ( (data[8].missing != -1) && (data[8].fvalue < (float)298752)) {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)237690)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)552)) {
            result[0] += -0.41628265;
          } else {
            result[0] += -0.17103384;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)197)) {
            result[0] += -0.34891665;
          } else {
            result[0] += 0.42154166;
          }
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)2058)) {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)100182)) {
            result[0] += 0.3135237;
          } else {
            result[0] += -0.49261427;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)117952)) {
            result[0] += 0.36800495;
          } else {
            result[0] += 0.056270074;
          }
        }
      }
    }
  } else {
    if ( (data[7].missing != -1) && (data[7].fvalue < (float)547712)) {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)45661)) {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)48946)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)47253)) {
            result[0] += 0.15922245;
          } else {
            result[0] += -0.19054948;
          }
        } else {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)343786)) {
            result[0] += -0.37296093;
          } else {
            result[0] += 0.10313491;
          }
        }
      } else {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)479164)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)879049)) {
            result[0] += -0.38585323;
          } else {
            result[0] += 0.3805061;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)139042)) {
            result[0] += -0.3937619;
          } else {
            result[0] += 0.3578938;
          }
        }
      }
    } else {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)323)) {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)124900)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)8229)) {
            result[0] += 0.35444745;
          } else {
            result[0] += -0.20431006;
          }
        } else {
          result[0] += -0.5795128;
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)319)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)94484)) {
            result[0] += 0.30274335;
          } else {
            result[0] += -0.72453374;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)503910)) {
            result[0] += 0.372095;
          } else {
            result[0] += -0.37470958;
          }
        }
      }
    }
  }
  if ( (data[22].missing != -1) && (data[22].fvalue < (float)814753)) {
    if ( (data[19].missing != -1) && (data[19].fvalue < (float)29865)) {
      if ( (data[19].missing != -1) && (data[19].fvalue < (float)418)) {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)95702)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)626)) {
            result[0] += -0.04704522;
          } else {
            result[0] += 0.14855762;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)18846)) {
            result[0] += -0.12122559;
          } else {
            result[0] += -0.42979965;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)57632)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)415)) {
            result[0] += 0.12394612;
          } else {
            result[0] += 0.35888478;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)48747)) {
            result[0] += 0.28972408;
          } else {
            result[0] += -0.07983421;
          }
        }
      }
    } else {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)211037)) {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)25658)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)676)) {
            result[0] += -0.14341693;
          } else {
            result[0] += 0.17077915;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)244910)) {
            result[0] += -0.32109895;
          } else {
            result[0] += 0.06177208;
          }
        }
      } else {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)250202)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)268)) {
            result[0] += -0.40453145;
          } else {
            result[0] += 0.17208807;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)432854)) {
            result[0] += -0.5125137;
          } else {
            result[0] += 0.24397127;
          }
        }
      }
    }
  } else {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)522)) {
      result[0] += -0.4586453;
    } else {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)212713)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)126)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)106615)) {
            result[0] += -0.5625138;
          } else {
            result[0] += 0.33776006;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)992579)) {
            result[0] += 0.23887502;
          } else {
            result[0] += 0.39561555;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)77696)) {
          result[0] += -0.7617165;
        } else {
          result[0] += 0.27196732;
        }
      }
    }
  }
  if ( (data[7].missing != -1) && (data[7].fvalue < (float)133)) {
    if ( (data[20].missing != -1) && (data[20].fvalue < (float)315)) {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)529)) {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)33788)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)158874)) {
            result[0] += -0.369018;
          } else {
            result[0] += -0.07933341;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)246154)) {
            result[0] += -0.102679566;
          } else {
            result[0] += 0.9794879;
          }
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)67560)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)150079)) {
            result[0] += -0.33574954;
          } else {
            result[0] += 0.029276941;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)351260)) {
            result[0] += 0.036476832;
          } else {
            result[0] += 1.1795459;
          }
        }
      }
    } else {
      if ( (data[14].missing != -1) && (data[14].fvalue < (float)152647)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)328390)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)50455)) {
            result[0] += 0.09136974;
          } else {
            result[0] += -0.056586243;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)138200)) {
            result[0] += 0.29734835;
          } else {
            result[0] += -0.58629936;
          }
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)476631)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)384)) {
            result[0] += 0.6094886;
          } else {
            result[0] += -0.48763928;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)392978)) {
            result[0] += 0.2789242;
          } else {
            result[0] += -0.6445992;
          }
        }
      }
    }
  } else {
    if ( (data[5].missing != -1) && (data[5].fvalue < (float)35548)) {
      if ( (data[5].missing != -1) && (data[5].fvalue < (float)635)) {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)78582)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)63)) {
            result[0] += -0.05622961;
          } else {
            result[0] += 0.16899873;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)25387)) {
            result[0] += -0.070318714;
          } else {
            result[0] += -0.32609007;
          }
        }
      } else {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)92611)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)193)) {
            result[0] += 0.22767203;
          } else {
            result[0] += 0.42820945;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)8640)) {
            result[0] += 0.28075793;
          } else {
            result[0] += -0.21632765;
          }
        }
      }
    } else {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)28312)) {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)706)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)114314)) {
            result[0] += 0.00895532;
          } else {
            result[0] += -0.35988078;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)89847)) {
            result[0] += 0.35934246;
          } else {
            result[0] += -0.0051973686;
          }
        }
      } else {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)572401)) {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)30422)) {
            result[0] += -0.07148668;
          } else {
            result[0] += -0.32316712;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)48202)) {
            result[0] += -0.45567232;
          } else {
            result[0] += 0.4162427;
          }
        }
      }
    }
  }
  if ( (data[23].missing != -1) && (data[23].fvalue < (float)909421)) {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)856828)) {
      if ( (data[24].missing != -1) && (data[24].fvalue < (float)57632)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)529)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)136372)) {
            result[0] += -0.08772215;
          } else {
            result[0] += 0.11571983;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)9791)) {
            result[0] += 0.20881024;
          } else {
            result[0] += 0.0505099;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)177737)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)25306)) {
            result[0] += 0.0035363608;
          } else {
            result[0] += -0.2702466;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)199744)) {
            result[0] += 0.12387311;
          } else {
            result[0] += -0.49370295;
          }
        }
      }
    } else {
      if ( (data[19].missing != -1) && (data[19].fvalue < (float)125298)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)147228)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)402)) {
            result[0] += -0.02554071;
          } else {
            result[0] += 0.46255705;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)240543)) {
            result[0] += -0.3393906;
          } else {
            result[0] += 0.17036031;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)431406)) {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)49160)) {
            result[0] += -0.6110062;
          } else {
            result[0] += 0.2200454;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)86542)) {
            result[0] += 0.3830851;
          } else {
            result[0] += 0.04857886;
          }
        }
      }
    }
  } else {
    if ( (data[24].missing != -1) && (data[24].fvalue < (float)9481)) {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)1159937)) {
        result[0] += -0.45836526;
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)6915)) {
          result[0] += -0.07424843;
        } else {
          result[0] += 0.31069207;
        }
      }
    } else {
      if ( (data[5].missing != -1) && (data[5].fvalue < (float)218493)) {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)190129)) {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)437883)) {
            result[0] += 0.41608652;
          } else {
            result[0] += 0.013873182;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)78073)) {
            result[0] += -0.13937253;
          } else {
            result[0] += 0.21489103;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)193301)) {
          result[0] += -0.25883225;
        } else {
          result[0] += 0.29045996;
        }
      }
    }
  }
  if ( (data[0].missing != -1) && (data[0].fvalue < (float)193)) {
    if ( (data[14].missing != -1) && (data[14].fvalue < (float)60904)) {
      if ( (data[14].missing != -1) && (data[14].fvalue < (float)196)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)64)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)220552)) {
            result[0] += -0.26123613;
          } else {
            result[0] += -0.73725635;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)217540)) {
            result[0] += -0.06573289;
          } else {
            result[0] += 0.14340684;
          }
        }
      } else {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)33883)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)127544)) {
            result[0] += 0.14683014;
          } else {
            result[0] += -0.15472375;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)579422)) {
            result[0] += -0.38047463;
          } else {
            result[0] += 0.3634338;
          }
        }
      }
    } else {
      if ( (data[8].missing != -1) && (data[8].fvalue < (float)197775)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)282509)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)154969)) {
            result[0] += -0.19872634;
          } else {
            result[0] += -0.43376657;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)88314)) {
            result[0] += 0.16121115;
          } else {
            result[0] += -0.33806276;
          }
        }
      } else {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)599)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)89263)) {
            result[0] += -0.72822684;
          } else {
            result[0] += 0.00032465399;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)60106)) {
            result[0] += 0.4425652;
          } else {
            result[0] += -0.0026172434;
          }
        }
      }
    }
  } else {
    if ( (data[0].missing != -1) && (data[0].fvalue < (float)48946)) {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)681)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)234)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)128)) {
            result[0] += -0.3445456;
          } else {
            result[0] += 0.01956593;
          }
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)114505)) {
            result[0] += 0.118908666;
          } else {
            result[0] += -0.58349866;
          }
        }
      } else {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)20274)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)221791)) {
            result[0] += 0.38338384;
          } else {
            result[0] += -0.039718978;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)204219)) {
            result[0] += 0.09454567;
          } else {
            result[0] += 0.33152595;
          }
        }
      }
    } else {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)19895)) {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)46972)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)286345)) {
            result[0] += 0.11432637;
          } else {
            result[0] += -0.53088504;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)1731)) {
            result[0] += -0.29811513;
          } else {
            result[0] += 0.055956185;
          }
        }
      } else {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)302106)) {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)124)) {
            result[0] += -0.055751026;
          } else {
            result[0] += -0.29299065;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)10319)) {
            result[0] += -0.34992728;
          } else {
            result[0] += 0.16662167;
          }
        }
      }
    }
  }
  if ( (data[12].missing != -1) && (data[12].fvalue < (float)421)) {
    if ( (data[17].missing != -1) && (data[17].fvalue < (float)60914)) {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)196)) {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)294083)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)222424)) {
            result[0] += -0.15078299;
          } else {
            result[0] += 0.15418084;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)70816)) {
            result[0] += 0.31239986;
          } else {
            result[0] += -0.20955186;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)129)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)300)) {
            result[0] += -0.3696351;
          } else {
            result[0] += 0.046299037;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)274961)) {
            result[0] += 0.16436079;
          } else {
            result[0] += -0.18140693;
          }
        }
      }
    } else {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)9062)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)39954)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)61)) {
            result[0] += -0.23237638;
          } else {
            result[0] += 0.065909944;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)691562)) {
            result[0] += -0.38550827;
          } else {
            result[0] += 0.44316193;
          }
        }
      } else {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)535979)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)76943)) {
            result[0] += -0.24166472;
          } else {
            result[0] += -0.47186166;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)131718)) {
            result[0] += 0.34903294;
          } else {
            result[0] += -0.6046817;
          }
        }
      }
    }
  } else {
    if ( (data[12].missing != -1) && (data[12].fvalue < (float)34330)) {
      if ( (data[17].missing != -1) && (data[17].fvalue < (float)55577)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)317)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)131)) {
            result[0] += -0.11650517;
          } else {
            result[0] += 0.22955282;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)3623)) {
            result[0] += 0.27382562;
          } else {
            result[0] += 0.52352995;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)75113)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)148080)) {
            result[0] += 0.012277792;
          } else {
            result[0] += 0.26493776;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)62837)) {
            result[0] += 0.031356096;
          } else {
            result[0] += -0.4180216;
          }
        }
      }
    } else {
      if ( (data[17].missing != -1) && (data[17].fvalue < (float)63)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)17146)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)54651)) {
            result[0] += 0.044815075;
          } else {
            result[0] += -0.22837281;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)166127)) {
            result[0] += -0.4374144;
          } else {
            result[0] += -0.07736537;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)95861)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)107646)) {
            result[0] += 0.09886158;
          } else {
            result[0] += -0.12352874;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)202035)) {
            result[0] += -0.27309567;
          } else {
            result[0] += 0.14804024;
          }
        }
      }
    }
  }
  if ( (data[4].missing != -1) && (data[4].fvalue < (float)216486)) {
    if ( (data[5].missing != -1) && (data[5].fvalue < (float)785224)) {
      if ( (data[0].missing != -1) && (data[0].fvalue < (float)193769)) {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)193)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)215018)) {
            result[0] += -0.017388681;
          } else {
            result[0] += -0.6609305;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)20059)) {
            result[0] += 0.13999727;
          } else {
            result[0] += 0.0012250885;
          }
        }
      } else {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)476425)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)66691)) {
            result[0] += -0.0816948;
          } else {
            result[0] += -0.3630621;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)314610)) {
            result[0] += 0.5183995;
          } else {
            result[0] += -0.44231406;
          }
        }
      }
    } else {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)238293)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)122105)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)273573)) {
            result[0] += -0.48163992;
          } else {
            result[0] += 0.027036378;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)61)) {
            result[0] += 0.47204772;
          } else {
            result[0] += 0.023504164;
          }
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)101281)) {
          result[0] += -0.3646734;
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)37123)) {
            result[0] += 0.7779307;
          } else {
            result[0] += 0.46476552;
          }
        }
      }
    }
  } else {
    if ( (data[3].missing != -1) && (data[3].fvalue < (float)629733)) {
      if ( (data[9].missing != -1) && (data[9].fvalue < (float)890202)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)13961)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)62)) {
            result[0] += -0.24087964;
          } else {
            result[0] += -0.5666097;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)7153)) {
            result[0] += -0.09248477;
          } else {
            result[0] += -0.33980998;
          }
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)7009)) {
          result[0] += -0.03732757;
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)137603)) {
            result[0] += 0.6789229;
          } else {
            result[0] += 0.21100178;
          }
        }
      }
    } else {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)225893)) {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)16356)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)44391)) {
            result[0] += -0.40651798;
          } else {
            result[0] += 0.23180854;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)442019)) {
            result[0] += 0.40271726;
          } else {
            result[0] += -0.34027588;
          }
        }
      } else {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)75751)) {
          result[0] += -0.83131516;
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)240754)) {
            result[0] += 0.26152295;
          } else {
            result[0] += 0.037579875;
          }
        }
      }
    }
  }
  if ( (data[18].missing != -1) && (data[18].fvalue < (float)733467)) {
    if ( (data[19].missing != -1) && (data[19].fvalue < (float)94033)) {
      if ( (data[14].missing != -1) && (data[14].fvalue < (float)112107)) {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)173237)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)12863)) {
            result[0] += 0.05732317;
          } else {
            result[0] += -0.06833711;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)603)) {
            result[0] += -0.074299134;
          } else {
            result[0] += 0.31186315;
          }
        }
      } else {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)63)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)38513)) {
            result[0] += -0.17971091;
          } else {
            result[0] += -0.43793398;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)14257)) {
            result[0] += 0.09840635;
          } else {
            result[0] += -0.17140295;
          }
        }
      }
    } else {
      if ( (data[8].missing != -1) && (data[8].fvalue < (float)195078)) {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)68068)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)163930)) {
            result[0] += -0.16947405;
          } else {
            result[0] += 0.1070685;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)513338)) {
            result[0] += -0.3870432;
          } else {
            result[0] += 0.37698182;
          }
        }
      } else {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)397)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)2433)) {
            result[0] += 0.35198885;
          } else {
            result[0] += -0.49716184;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)85209)) {
            result[0] += 0.34494826;
          } else {
            result[0] += 0.008838673;
          }
        }
      }
    }
  } else {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)656)) {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)456649)) {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)166455)) {
          result[0] += -0.61562425;
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)86954)) {
            result[0] += 0.14024097;
          } else {
            result[0] += -0.32234725;
          }
        }
      } else {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)137603)) {
          result[0] += -0.25259644;
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)7602)) {
            result[0] += 0.03890536;
          } else {
            result[0] += 0.33704564;
          }
        }
      }
    } else {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)125)) {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)89395)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)18511)) {
            result[0] += -0.025281487;
          } else {
            result[0] += -0.7432101;
          }
        } else {
          result[0] += 0.5225317;
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)61)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)61)) {
            result[0] += 0.22316042;
          } else {
            result[0] += -0.49455795;
          }
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)252153)) {
            result[0] += 0.38189897;
          } else {
            result[0] += -0.018927926;
          }
        }
      }
    }
  }
  if ( (data[4].missing != -1) && (data[4].fvalue < (float)140493)) {
    if ( (data[3].missing != -1) && (data[3].fvalue < (float)846)) {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)253776)) {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)129)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)706)) {
            result[0] += -0.26186857;
          } else {
            result[0] += -0.035900764;
          }
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)142324)) {
            result[0] += 0.051491667;
          } else {
            result[0] += -0.34332028;
          }
        }
      } else {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)273766)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)452648)) {
            result[0] += -0.82568115;
          } else {
            result[0] += -0.1776233;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)225975)) {
            result[0] += -0.31764668;
          } else {
            result[0] += 0.5953638;
          }
        }
      }
    } else {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)21425)) {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)191006)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)130)) {
            result[0] += 0.07208108;
          } else {
            result[0] += 0.2685084;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)944)) {
            result[0] += -0.5019554;
          } else {
            result[0] += 0.107112005;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)44394)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)2758)) {
            result[0] += 0.14649712;
          } else {
            result[0] += -0.034933407;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)25422)) {
            result[0] += 0.1137295;
          } else {
            result[0] += -0.14819096;
          }
        }
      }
    }
  } else {
    if ( (data[15].missing != -1) && (data[15].fvalue < (float)33464)) {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)127)) {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)208874)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)102337)) {
            result[0] += -0.27429757;
          } else {
            result[0] += 0.14037596;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)642320)) {
            result[0] += -0.5573794;
          } else {
            result[0] += -0.09683921;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)1326)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)1576)) {
            result[0] += -0.4886342;
          } else {
            result[0] += 0.11306182;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)25306)) {
            result[0] += 0.5747691;
          } else {
            result[0] += 0.108979166;
          }
        }
      }
    } else {
      if ( (data[8].missing != -1) && (data[8].fvalue < (float)709386)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)341687)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)711838)) {
            result[0] += -0.41823658;
          } else {
            result[0] += 0.13098203;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)134880)) {
            result[0] += -0.3093513;
          } else {
            result[0] += 0.33972302;
          }
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)849)) {
          result[0] += -0.29889795;
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)692)) {
            result[0] += -0.048612718;
          } else {
            result[0] += 0.41768643;
          }
        }
      }
    }
  }
  if ( (data[10].missing != -1) && (data[10].fvalue < (float)35064)) {
    if ( (data[10].missing != -1) && (data[10].fvalue < (float)563)) {
      if ( (data[7].missing != -1) && (data[7].fvalue < (float)92236)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)106132)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)446)) {
            result[0] += -0.059072297;
          } else {
            result[0] += 0.13510697;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)278748)) {
            result[0] += -0.18646003;
          } else {
            result[0] += 0.19178478;
          }
        }
      } else {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)105074)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)2978)) {
            result[0] += -0.22260956;
          } else {
            result[0] += 0.08425247;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)373376)) {
            result[0] += -0.5124768;
          } else {
            result[0] += 0.2821059;
          }
        }
      }
    } else {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)95433)) {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)280080)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)159848)) {
            result[0] += 0.17029296;
          } else {
            result[0] += 0.38491854;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)136372)) {
            result[0] += -0.41745374;
          } else {
            result[0] += 0.34080088;
          }
        }
      } else {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)187721)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)403367)) {
            result[0] += -0.12517917;
          } else {
            result[0] += 0.3106633;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)361937)) {
            result[0] += 0.384945;
          } else {
            result[0] += 0.041404407;
          }
        }
      }
    }
  } else {
    if ( (data[6].missing != -1) && (data[6].fvalue < (float)467)) {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)68140)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)139276)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)113167)) {
            result[0] += -0.08679242;
          } else {
            result[0] += -0.33028606;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)46515)) {
            result[0] += 0.20664068;
          } else {
            result[0] += -0.29208162;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)718657)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)5101)) {
            result[0] += -0.5316237;
          } else {
            result[0] += -0.3125988;
          }
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)62)) {
            result[0] += 0.36992213;
          } else {
            result[0] += 0.08523325;
          }
        }
      }
    } else {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)28157)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)525)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)64222)) {
            result[0] += 0.08241349;
          } else {
            result[0] += -0.2429241;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)592770)) {
            result[0] += 0.23485938;
          } else {
            result[0] += 0.6487533;
          }
        }
      } else {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)624476)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)221656)) {
            result[0] += -0.21741365;
          } else {
            result[0] += 0.06816964;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)8640)) {
            result[0] += -0.592894;
          } else {
            result[0] += 0.4260442;
          }
        }
      }
    }
  }
  if ( (data[9].missing != -1) && (data[9].fvalue < (float)70235)) {
    if ( (data[22].missing != -1) && (data[22].fvalue < (float)191)) {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)569)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)97141)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)129403)) {
            result[0] += -0.2879611;
          } else {
            result[0] += -0.00036344436;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)55887)) {
            result[0] += -0.025296474;
          } else {
            result[0] += 0.4850761;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)96243)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)175977)) {
            result[0] += 0.07384111;
          } else {
            result[0] += 0.38856465;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)367906)) {
            result[0] += -0.2508199;
          } else {
            result[0] += 0.20664985;
          }
        }
      }
    } else {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)50455)) {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)202)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)664)) {
            result[0] += -0.10821457;
          } else {
            result[0] += 0.11922147;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)237422)) {
            result[0] += 0.24230497;
          } else {
            result[0] += -0.1821981;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)324266)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)26794)) {
            result[0] += 0.118789375;
          } else {
            result[0] += -0.08455505;
          }
        } else {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)228995)) {
            result[0] += 0.18435945;
          } else {
            result[0] += -0.7332028;
          }
        }
      }
    }
  } else {
    if ( (data[3].missing != -1) && (data[3].fvalue < (float)156911)) {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)83716)) {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)79383)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)83257)) {
            result[0] += 0.08165507;
          } else {
            result[0] += -0.20939127;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)816)) {
            result[0] += -0.17525597;
          } else {
            result[0] += -0.4586528;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)468884)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)83147)) {
            result[0] += -0.08268427;
          } else {
            result[0] += -0.35588324;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)557)) {
            result[0] += -0.4392836;
          } else {
            result[0] += 0.19621876;
          }
        }
      }
    } else {
      if ( (data[14].missing != -1) && (data[14].fvalue < (float)598)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)507221)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)74885)) {
            result[0] += -0.15902491;
          } else {
            result[0] += -0.60616964;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)197)) {
            result[0] += -0.43933636;
          } else {
            result[0] += 0.5841854;
          }
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)2180)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)19571)) {
            result[0] += -0.019743664;
          } else {
            result[0] += -0.485382;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)260466)) {
            result[0] += 0.30815825;
          } else {
            result[0] += -0.28692207;
          }
        }
      }
    }
  }
  if ( (data[20].missing != -1) && (data[20].fvalue < (float)32652)) {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)522)) {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)294)) {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)137553)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)93372)) {
            result[0] += -0.24427131;
          } else {
            result[0] += -0.028283475;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)105972)) {
            result[0] += 0.23040819;
          } else {
            result[0] += -0.36923116;
          }
        }
      } else {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)16982)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)490)) {
            result[0] += 0.058136124;
          } else {
            result[0] += 0.35337695;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)567)) {
            result[0] += -0.20471536;
          } else {
            result[0] += -0.003973923;
          }
        }
      }
    } else {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)119952)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)1136)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)191)) {
            result[0] += -0.052252073;
          } else {
            result[0] += 0.17777425;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)123853)) {
            result[0] += 0.32291988;
          } else {
            result[0] += 0.06359827;
          }
        }
      } else {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)29362)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)234777)) {
            result[0] += 0.14349014;
          } else {
            result[0] += -0.54128534;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)855)) {
            result[0] += 0.52510834;
          } else {
            result[0] += -0.2606195;
          }
        }
      }
    }
  } else {
    if ( (data[14].missing != -1) && (data[14].fvalue < (float)39954)) {
      if ( (data[9].missing != -1) && (data[9].fvalue < (float)107919)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)200078)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)462239)) {
            result[0] += -0.024437232;
          } else {
            result[0] += 0.27469063;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)90997)) {
            result[0] += 0.2625314;
          } else {
            result[0] += -0.1389799;
          }
        }
      } else {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)31088)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)62881)) {
            result[0] += 0.24911356;
          } else {
            result[0] += -0.20126854;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)305559)) {
            result[0] += -0.49673063;
          } else {
            result[0] += -0.08227151;
          }
        }
      }
    } else {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)413438)) {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)150453)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)81901)) {
            result[0] += -0.17006697;
          } else {
            result[0] += -0.39146706;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)90579)) {
            result[0] += 0.09959292;
          } else {
            result[0] += -0.20403711;
          }
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)191213)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)7220)) {
            result[0] += 0.14289726;
          } else {
            result[0] += -0.40909466;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)67)) {
            result[0] += -0.40596366;
          } else {
            result[0] += 0.33724496;
          }
        }
      }
    }
  }
  if ( (data[20].missing != -1) && (data[20].fvalue < (float)885135)) {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)25306)) {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)2510)) {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)282892)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)148080)) {
            result[0] += -0.014609692;
          } else {
            result[0] += 0.17694442;
          }
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)9640)) {
            result[0] += -0.1637172;
          } else {
            result[0] += -0.5622474;
          }
        }
      } else {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)329)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)109111)) {
            result[0] += 0.11961501;
          } else {
            result[0] += -0.4690581;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)45709)) {
            result[0] += 0.34561118;
          } else {
            result[0] += 0.13187693;
          }
        }
      }
    } else {
      if ( (data[19].missing != -1) && (data[19].fvalue < (float)63893)) {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)327)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)118450)) {
            result[0] += -0.04138224;
          } else {
            result[0] += 0.11851525;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)35258)) {
            result[0] += 0.0076544033;
          } else {
            result[0] += -0.1780896;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)557222)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)456649)) {
            result[0] += -0.18444617;
          } else {
            result[0] += 0.18006073;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)4686)) {
            result[0] += -0.3011162;
          } else {
            result[0] += 0.2761474;
          }
        }
      }
    }
  } else {
    if ( (data[17].missing != -1) && (data[17].fvalue < (float)62)) {
      result[0] += -0.77749056;
    } else {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)29970)) {
        result[0] += -0.3891581;
      } else {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)144073)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)497537)) {
            result[0] += 0.59800565;
          } else {
            result[0] += 0.3727215;
          }
        } else {
          result[0] += -0.13928947;
        }
      }
    }
  }
  if ( (data[23].missing != -1) && (data[23].fvalue < (float)529)) {
    if ( (data[15].missing != -1) && (data[15].fvalue < (float)2506)) {
      if ( (data[24].missing != -1) && (data[24].fvalue < (float)61)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)257096)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)276)) {
            result[0] += -0.312671;
          } else {
            result[0] += -0.1049323;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)225)) {
            result[0] += 0.44867274;
          } else {
            result[0] += -0.3359197;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)471970)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)2293)) {
            result[0] += -0.09409318;
          } else {
            result[0] += 0.105043724;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)61)) {
            result[0] += 1.6882515;
          } else {
            result[0] += -0.09852525;
          }
        }
      }
    } else {
      if ( (data[9].missing != -1) && (data[9].fvalue < (float)96907)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)297045)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)498)) {
            result[0] += -0.032723937;
          } else {
            result[0] += 0.1411919;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)43068)) {
            result[0] += -0.39250377;
          } else {
            result[0] += 0.06662669;
          }
        }
      } else {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)162784)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)243969)) {
            result[0] += -0.32478628;
          } else {
            result[0] += 0.011483396;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)9189)) {
            result[0] += 0.3480598;
          } else {
            result[0] += -0.18030314;
          }
        }
      }
    }
  } else {
    if ( (data[23].missing != -1) && (data[23].fvalue < (float)18890)) {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)322)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)134)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)284912)) {
            result[0] += -0.32786402;
          } else {
            result[0] += 0.4235807;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)189)) {
            result[0] += 1.5751652;
          } else {
            result[0] += 0.06777746;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)402)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)557)) {
            result[0] += -0.050939865;
          } else {
            result[0] += 0.19677739;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)127691)) {
            result[0] += 0.33446774;
          } else {
            result[0] += -0.049779635;
          }
        }
      }
    } else {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)42744)) {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)128902)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)178619)) {
            result[0] += 0.10323238;
          } else {
            result[0] += -0.29102442;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)100644)) {
            result[0] += -0.018575614;
          } else {
            result[0] += -0.42760798;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)196362)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)106132)) {
            result[0] += -0.07324129;
          } else {
            result[0] += -0.2989554;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)438)) {
            result[0] += -0.2353094;
          } else {
            result[0] += 0.108338095;
          }
        }
      }
    }
  }
  if ( (data[21].missing != -1) && (data[21].fvalue < (float)987697)) {
    if ( (data[23].missing != -1) && (data[23].fvalue < (float)909421)) {
      if ( (data[11].missing != -1) && (data[11].fvalue < (float)615)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)266067)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)277786)) {
            result[0] += -0.021161282;
          } else {
            result[0] += -0.34468636;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)2077)) {
            result[0] += -0.75558656;
          } else {
            result[0] += -0.185393;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)18447)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)439150)) {
            result[0] += 0.10302361;
          } else {
            result[0] += -0.47040668;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)136862)) {
            result[0] += -0.0630239;
          } else {
            result[0] += 0.05924744;
          }
        }
      }
    } else {
      if ( (data[24].missing != -1) && (data[24].fvalue < (float)11975)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)303)) {
          result[0] += -0.7715959;
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)322)) {
            result[0] += 0.26891878;
          } else {
            result[0] += -0.103643894;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)1384)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)61)) {
            result[0] += -0.62722945;
          } else {
            result[0] += 0.29735082;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)249566)) {
            result[0] += 0.40083465;
          } else {
            result[0] += 0.008566014;
          }
        }
      }
    }
  } else {
    if ( (data[14].missing != -1) && (data[14].fvalue < (float)232154)) {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)60068)) {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)81901)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)247330)) {
            result[0] += 0.4627343;
          } else {
            result[0] += 0.14003505;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)77842)) {
            result[0] += -0.21110228;
          } else {
            result[0] += 0.25160718;
          }
        }
      } else {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)47196)) {
          result[0] += 0.30818373;
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)21319)) {
            result[0] += 0.03305194;
          } else {
            result[0] += -0.42773315;
          }
        }
      }
    } else {
      result[0] += -0.22842343;
    }
  }
  if ( (data[22].missing != -1) && (data[22].fvalue < (float)1100503)) {
    if ( (data[1].missing != -1) && (data[1].fvalue < (float)301165)) {
      if ( (data[1].missing != -1) && (data[1].fvalue < (float)596)) {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)261261)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)291392)) {
            result[0] += -0.0016570495;
          } else {
            result[0] += -0.54657984;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)94863)) {
            result[0] += -0.66529596;
          } else {
            result[0] += -0.06222671;
          }
        }
      } else {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)33883)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)285)) {
            result[0] += 0.032040328;
          } else {
            result[0] += 0.16942935;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)17146)) {
            result[0] += 0.0612301;
          } else {
            result[0] += -0.08757415;
          }
        }
      }
    } else {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)123670)) {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)1482)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)70459)) {
            result[0] += -0.35452494;
          } else {
            result[0] += 0.41773134;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)988220)) {
            result[0] += -0.6947021;
          } else {
            result[0] += -0.101329006;
          }
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)199166)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)63989)) {
            result[0] += -0.068872;
          } else {
            result[0] += -0.59054667;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)61)) {
            result[0] += -0.3103225;
          } else {
            result[0] += 0.2667065;
          }
        }
      }
    }
  } else {
    if ( (data[16].missing != -1) && (data[16].fvalue < (float)1022)) {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)11842)) {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)36635)) {
          result[0] += 0.33650962;
        } else {
          result[0] += 0.020306004;
        }
      } else {
        result[0] += -0.7846388;
      }
    } else {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)452550)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)4716)) {
          result[0] += -0.035459843;
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)150996)) {
            result[0] += 0.4061072;
          } else {
            result[0] += 0.12815768;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)625770)) {
          result[0] += -0.39858875;
        } else {
          result[0] += 0.32530165;
        }
      }
    }
  }
  if ( (data[17].missing != -1) && (data[17].fvalue < (float)620)) {
    if ( (data[6].missing != -1) && (data[6].fvalue < (float)131)) {
      if ( (data[0].missing != -1) && (data[0].fvalue < (float)107512)) {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)314)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)35823)) {
            result[0] += -0.38686788;
          } else {
            result[0] += -0.14851233;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)52110)) {
            result[0] += -0.12554698;
          } else {
            result[0] += 0.2051479;
          }
        }
      } else {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)110013)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)62561)) {
            result[0] += 0.03605927;
          } else {
            result[0] += 1.309948;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)416)) {
            result[0] += -0.12736312;
          } else {
            result[0] += 0.22946782;
          }
        }
      }
    } else {
      if ( (data[14].missing != -1) && (data[14].fvalue < (float)296398)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)352714)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)159459)) {
            result[0] += -0.0038365491;
          } else {
            result[0] += 0.23611097;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)277786)) {
            result[0] += -0.5337241;
          } else {
            result[0] += 0.35660207;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)64)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)404783)) {
            result[0] += -0.46203613;
          } else {
            result[0] += -0.77816;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)13272)) {
            result[0] += 0.027611697;
          } else {
            result[0] += -0.47618452;
          }
        }
      }
    }
  } else {
    if ( (data[16].missing != -1) && (data[16].fvalue < (float)36501)) {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)65)) {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)119216)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)721)) {
            result[0] += -0.10923512;
          } else {
            result[0] += 0.15760288;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)168443)) {
            result[0] += -0.42422614;
          } else {
            result[0] += 0.06926011;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)143078)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)349425)) {
            result[0] += 0.23430093;
          } else {
            result[0] += -0.21249382;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)73823)) {
            result[0] += 0.08265209;
          } else {
            result[0] += -0.3703108;
          }
        }
      }
    } else {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)144953)) {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)43068)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)155715)) {
            result[0] += 0.0060480516;
          } else {
            result[0] += 0.31607124;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)22939)) {
            result[0] += -0.033398777;
          } else {
            result[0] += -0.17578118;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)95238)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)265244)) {
            result[0] += 0.19865845;
          } else {
            result[0] += -0.24953063;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)323290)) {
            result[0] += -0.13736662;
          } else {
            result[0] += 0.24286032;
          }
        }
      }
    }
  }
  if ( (data[6].missing != -1) && (data[6].fvalue < (float)44330)) {
    if ( (data[6].missing != -1) && (data[6].fvalue < (float)544)) {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)170570)) {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)276112)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)134732)) {
            result[0] += 0.00551347;
          } else {
            result[0] += -0.1525623;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)498)) {
            result[0] += -0.41755113;
          } else {
            result[0] += 0.18291548;
          }
        }
      } else {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)3977)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)138484)) {
            result[0] += -0.7328462;
          } else {
            result[0] += -0.21441896;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)471)) {
            result[0] += 0.05218458;
          } else {
            result[0] += -0.3327825;
          }
        }
      }
    } else {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)748)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)212)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)62492)) {
            result[0] += -0.0075923293;
          } else {
            result[0] += 0.25571144;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)371270)) {
            result[0] += -0.6075935;
          } else {
            result[0] += 0.21738124;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)553)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)136342)) {
            result[0] += 0.1462131;
          } else {
            result[0] += -0.1827465;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)9743)) {
            result[0] += 0.36601713;
          } else {
            result[0] += 0.16736451;
          }
        }
      }
    }
  } else {
    if ( (data[4].missing != -1) && (data[4].fvalue < (float)9832)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)859)) {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)3889)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)173838)) {
            result[0] += 0.078790404;
          } else {
            result[0] += -0.40480945;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)575825)) {
            result[0] += -0.1791016;
          } else {
            result[0] += 0.42307195;
          }
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)65173)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)116067)) {
            result[0] += -0.42587286;
          } else {
            result[0] += 0.19016643;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)89219)) {
            result[0] += 0.48159343;
          } else {
            result[0] += 0.16388175;
          }
        }
      }
    } else {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)893)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)454241)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)1490932)) {
            result[0] += -0.43662128;
          } else {
            result[0] += 0.44768786;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)81901)) {
            result[0] += -0.35097063;
          } else {
            result[0] += 0.33368394;
          }
        }
      } else {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)237690)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)244910)) {
            result[0] += -0.1957213;
          } else {
            result[0] += 0.08055619;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)337812)) {
            result[0] += 0.18079299;
          } else {
            result[0] += -0.33534387;
          }
        }
      }
    }
  }
  if ( (data[5].missing != -1) && (data[5].fvalue < (float)25658)) {
    if ( (data[5].missing != -1) && (data[5].fvalue < (float)528)) {
      if ( (data[8].missing != -1) && (data[8].fvalue < (float)101511)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)4268)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)62)) {
            result[0] += -0.218238;
          } else {
            result[0] += -0.0043309494;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)27198)) {
            result[0] += 0.13142832;
          } else {
            result[0] += -0.023423146;
          }
        }
      } else {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)164974)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)352030)) {
            result[0] += -0.2279177;
          } else {
            result[0] += -0.6584127;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)83937)) {
            result[0] += -0.19535272;
          } else {
            result[0] += 0.2822591;
          }
        }
      }
    } else {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)265244)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)180247)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)32074)) {
            result[0] += 0.18091968;
          } else {
            result[0] += 0.024107467;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)111800)) {
            result[0] += 0.38879463;
          } else {
            result[0] += 0.019407824;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)45526)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)16311)) {
            result[0] += -0.63189775;
          } else {
            result[0] += -0.21230076;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)195)) {
            result[0] += -0.6403143;
          } else {
            result[0] += 0.16388862;
          }
        }
      }
    }
  } else {
    if ( (data[10].missing != -1) && (data[10].fvalue < (float)879049)) {
      if ( (data[11].missing != -1) && (data[11].fvalue < (float)26385)) {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)77842)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)71343)) {
            result[0] += 0.022908002;
          } else {
            result[0] += 0.58007;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)99133)) {
            result[0] += -0.08839088;
          } else {
            result[0] += -0.38777176;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)552)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)17633)) {
            result[0] += -0.17697015;
          } else {
            result[0] += -0.4937581;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)821)) {
            result[0] += -0.18425491;
          } else {
            result[0] += 0.00079251424;
          }
        }
      }
    } else {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)656)) {
        result[0] += -0.54105496;
      } else {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)182803)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)9743)) {
            result[0] += -0.1495615;
          } else {
            result[0] += 0.5234762;
          }
        } else {
          result[0] += -0.19403845;
        }
      }
    }
  }
  if ( (data[20].missing != -1) && (data[20].fvalue < (float)597342)) {
    if ( (data[18].missing != -1) && (data[18].fvalue < (float)224)) {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)301276)) {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)57642)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)235222)) {
            result[0] += -0.053416487;
          } else {
            result[0] += -0.5030507;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)115782)) {
            result[0] += 0.13526908;
          } else {
            result[0] += -0.24316342;
          }
        }
      } else {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)215)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)54651)) {
            result[0] += -0.7941746;
          } else {
            result[0] += -0.14318053;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)22540)) {
            result[0] += -0.024767054;
          } else {
            result[0] += -0.65615535;
          }
        }
      }
    } else {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)31207)) {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)273934)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)42617)) {
            result[0] += 0.13257375;
          } else {
            result[0] += 0.010849193;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)235132)) {
            result[0] += -0.3585792;
          } else {
            result[0] += 0.14637236;
          }
        }
      } else {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)55445)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)124634)) {
            result[0] += 0.03592251;
          } else {
            result[0] += -0.24287622;
          }
        } else {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)283270)) {
            result[0] += -0.14278002;
          } else {
            result[0] += 0.2375007;
          }
        }
      }
    }
  } else {
    if ( (data[18].missing != -1) && (data[18].fvalue < (float)157890)) {
      if ( (data[0].missing != -1) && (data[0].fvalue < (float)86275)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)222145)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)61)) {
            result[0] += 0.5128987;
          } else {
            result[0] += 0.22885731;
          }
        } else {
          result[0] += -0.4280579;
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)187953)) {
          result[0] += -0.6046047;
        } else {
          result[0] += 0.3508568;
        }
      }
    } else {
      if ( (data[17].missing != -1) && (data[17].fvalue < (float)718657)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)224313)) {
          result[0] += -0.6358648;
        } else {
          result[0] += 0.035073217;
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)1292)) {
          result[0] += -0.092759326;
        } else {
          result[0] += 0.33453605;
        }
      }
    }
  }
  if ( (data[16].missing != -1) && (data[16].fvalue < (float)61)) {
    if ( (data[24].missing != -1) && (data[24].fvalue < (float)159892)) {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)158038)) {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)855)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)614)) {
            result[0] += -0.20603676;
          } else {
            result[0] += -0.034946673;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)79364)) {
            result[0] += 0.047972135;
          } else {
            result[0] += -0.18560585;
          }
        }
      } else {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)239486)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)120792)) {
            result[0] += 0.23927946;
          } else {
            result[0] += -0.22042398;
          }
        } else {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)334070)) {
            result[0] += -0.5759924;
          } else {
            result[0] += 0.3349156;
          }
        }
      }
    } else {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)10659)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)104935)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)52732)) {
            result[0] += -0.3453892;
          } else {
            result[0] += 0.36141673;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)7134)) {
            result[0] += 0.19746302;
          } else {
            result[0] += -0.5387524;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)882427)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)747401)) {
            result[0] += -0.5468951;
          } else {
            result[0] += 0.026867157;
          }
        } else {
          result[0] += 0.38463634;
        }
      }
    }
  } else {
    if ( (data[16].missing != -1) && (data[16].fvalue < (float)14327)) {
      if ( (data[1].missing != -1) && (data[1].fvalue < (float)301165)) {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)299677)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)1090)) {
            result[0] += 0.03363117;
          } else {
            result[0] += 0.1634803;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)2064)) {
            result[0] += -0.40300712;
          } else {
            result[0] += 0.05306637;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)890202)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)693580)) {
            result[0] += -0.5975339;
          } else {
            result[0] += 0.06889715;
          }
        } else {
          result[0] += 0.16945295;
        }
      }
    } else {
      if ( (data[18].missing != -1) && (data[18].fvalue < (float)556445)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)33740)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)875)) {
            result[0] += -0.037960146;
          } else {
            result[0] += 0.11023742;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)63233)) {
            result[0] += 0.041734274;
          } else {
            result[0] += -0.15543759;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)1667)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)429)) {
            result[0] += -0.41214785;
          } else {
            result[0] += 0.11392546;
          }
        } else {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)135217)) {
            result[0] += 0.33933088;
          } else {
            result[0] += 0.07776737;
          }
        }
      }
    }
  }
  if ( (data[5].missing != -1) && (data[5].fvalue < (float)61)) {
    if ( (data[7].missing != -1) && (data[7].fvalue < (float)144611)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)108554)) {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)107499)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)500)) {
            result[0] += -0.084103905;
          } else {
            result[0] += 0.028136952;
          }
        } else {
          result[0] += 1.2107695;
        }
      } else {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)690647)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)347087)) {
            result[0] += -0.35499635;
          } else {
            result[0] += 0.17027414;
          }
        } else {
          result[0] += 0.766715;
        }
      }
    } else {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)300)) {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)1759)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)94515)) {
            result[0] += -1.0918818;
          } else {
            result[0] += -0.3136222;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)34957)) {
            result[0] += -0.206333;
          } else {
            result[0] += -0.65558416;
          }
        }
      } else {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)21824)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)532499)) {
            result[0] += -0.570822;
          } else {
            result[0] += 0.44869095;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)78730)) {
            result[0] += 0.19499068;
          } else {
            result[0] += -0.2825368;
          }
        }
      }
    }
  } else {
    if ( (data[5].missing != -1) && (data[5].fvalue < (float)20274)) {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)351214)) {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)249365)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)193237)) {
            result[0] += 0.08007866;
          } else {
            result[0] += 0.21970169;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)412089)) {
            result[0] += -0.32305315;
          } else {
            result[0] += 0.28148937;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)173522)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)339281)) {
            result[0] += -0.6682014;
          } else {
            result[0] += 0.0550672;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)361)) {
            result[0] += -0.45506343;
          } else {
            result[0] += 0.2842188;
          }
        }
      }
    } else {
      if ( (data[7].missing != -1) && (data[7].fvalue < (float)468884)) {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)32817)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)10405)) {
            result[0] += 0.07408248;
          } else {
            result[0] += -0.09664377;
          }
        } else {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)666602)) {
            result[0] += -0.10460738;
          } else {
            result[0] += 0.37050968;
          }
        }
      } else {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)164262)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)172008)) {
            result[0] += 0.2702185;
          } else {
            result[0] += -0.09822846;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)288822)) {
            result[0] += -0.4464762;
          } else {
            result[0] += 0.30234143;
          }
        }
      }
    }
  }
  if ( (data[16].missing != -1) && (data[16].fvalue < (float)61)) {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)136406)) {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)116647)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)78890)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)894)) {
            result[0] += -0.24062108;
          } else {
            result[0] += -0.0052402914;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)271436)) {
            result[0] += 0.11955642;
          } else {
            result[0] += -0.6793396;
          }
        }
      } else {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)214633)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)252153)) {
            result[0] += 0.18470117;
          } else {
            result[0] += -0.45897385;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)284596)) {
            result[0] += -0.630814;
          } else {
            result[0] += -0.13077025;
          }
        }
      }
    } else {
      if ( (data[11].missing != -1) && (data[11].fvalue < (float)189)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)218985)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)163586)) {
            result[0] += -0.3540261;
          } else {
            result[0] += -0.01707655;
          }
        } else {
          result[0] += -0.7629565;
        }
      } else {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)78663)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)62)) {
            result[0] += -0.3258317;
          } else {
            result[0] += -0.67212737;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)112970)) {
            result[0] += -0.018158242;
          } else {
            result[0] += -0.50986695;
          }
        }
      }
    }
  } else {
    if ( (data[4].missing != -1) && (data[4].fvalue < (float)16982)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)15066)) {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)391925)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)302)) {
            result[0] += 0.0058194;
          } else {
            result[0] += 0.07398829;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)18497)) {
            result[0] += -0.7907905;
          } else {
            result[0] += 0.029236097;
          }
        }
      } else {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)65)) {
          result[0] += -0.4079511;
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)254252)) {
            result[0] += 0.6049897;
          } else {
            result[0] += -0.36506054;
          }
        }
      }
    } else {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)46839)) {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)123853)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)18744)) {
            result[0] += 0.1592965;
          } else {
            result[0] += -0.03554168;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)125)) {
            result[0] += 0.7611729;
          } else {
            result[0] += -0.28818226;
          }
        }
      } else {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)501842)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)48120)) {
            result[0] += -0.063584335;
          } else {
            result[0] += -0.25685927;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)555)) {
            result[0] += -0.42146918;
          } else {
            result[0] += 0.23355661;
          }
        }
      }
    }
  }
  if ( (data[12].missing != -1) && (data[12].fvalue < (float)179982)) {
    if ( (data[9].missing != -1) && (data[9].fvalue < (float)20029)) {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)912)) {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)315)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)2496)) {
            result[0] += -0.20803152;
          } else {
            result[0] += -0.007093075;
          }
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)240754)) {
            result[0] += 0.0055239713;
          } else {
            result[0] += 0.35839888;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)15006)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)154969)) {
            result[0] += 0.09596721;
          } else {
            result[0] += -0.18654339;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)235)) {
            result[0] += -0.11986976;
          } else {
            result[0] += 0.051692363;
          }
        }
      }
    } else {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)166304)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)9791)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)748)) {
            result[0] += -0.1616671;
          } else {
            result[0] += 0.06328722;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)34330)) {
            result[0] += -0.117292486;
          } else {
            result[0] += -0.33187208;
          }
        }
      } else {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)2402)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)336637)) {
            result[0] += -0.5017281;
          } else {
            result[0] += 0.40214595;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)97773)) {
            result[0] += 0.20237653;
          } else {
            result[0] += -0.08630309;
          }
        }
      }
    }
  } else {
    if ( (data[12].missing != -1) && (data[12].fvalue < (float)247668)) {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)47307)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)161400)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)208238)) {
            result[0] += 0.5382865;
          } else {
            result[0] += -0.1525352;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)55376)) {
            result[0] += 0.10596102;
          } else {
            result[0] += -0.47289968;
          }
        }
      } else {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)65)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)554)) {
            result[0] += 0.42003784;
          } else {
            result[0] += -0.078367576;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)381251)) {
            result[0] += -0.19067682;
          } else {
            result[0] += 0.23658119;
          }
        }
      }
    } else {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)30345)) {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)312845)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)202339)) {
            result[0] += -0.020675171;
          } else {
            result[0] += -0.4904166;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)356962)) {
            result[0] += -0.7078241;
          } else {
            result[0] += 0.1434127;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)106904)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)136602)) {
            result[0] += 0.17405705;
          } else {
            result[0] += -0.25195053;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)555)) {
            result[0] += -0.52570695;
          } else {
            result[0] += -0.010885249;
          }
        }
      }
    }
  }
  if ( (data[22].missing != -1) && (data[22].fvalue < (float)1100503)) {
    if ( (data[2].missing != -1) && (data[2].fvalue < (float)833588)) {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)1282683)) {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)859444)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)108656)) {
            result[0] += 0.0033529743;
          } else {
            result[0] += -0.061462816;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)11848)) {
            result[0] += -0.63446426;
          } else {
            result[0] += 0.4577919;
          }
        }
      } else {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)1404)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)88052)) {
            result[0] += -0.38435712;
          } else {
            result[0] += 0.15599597;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)61023)) {
            result[0] += 0.5170088;
          } else {
            result[0] += -0.044771772;
          }
        }
      }
    } else {
      if ( (data[1].missing != -1) && (data[1].fvalue < (float)62246)) {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)55284)) {
          result[0] += -0.7194396;
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)118944)) {
            result[0] += -0.1293829;
          } else {
            result[0] += 0.29704604;
          }
        }
      } else {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)434)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)67655)) {
            result[0] += -0.2925258;
          } else {
            result[0] += 0.044648807;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)395684)) {
            result[0] += 0.39070296;
          } else {
            result[0] += -0.055065986;
          }
        }
      }
    }
  } else {
    if ( (data[1].missing != -1) && (data[1].fvalue < (float)79383)) {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)108066)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)15040)) {
          result[0] += -0.028608544;
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)622776)) {
            result[0] += 0.4033382;
          } else {
            result[0] += 0.064019956;
          }
        }
      } else {
        result[0] += -0.082220845;
      }
    } else {
      if ( (data[19].missing != -1) && (data[19].fvalue < (float)117919)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)9189)) {
          result[0] += 0.32912132;
        } else {
          result[0] += -0.17114364;
        }
      } else {
        result[0] += -0.5363989;
      }
    }
  }
  if ( (data[7].missing != -1) && (data[7].fvalue < (float)133)) {
    if ( (data[7].missing != -1) && (data[7].fvalue < (float)129)) {
      if ( (data[11].missing != -1) && (data[11].fvalue < (float)121285)) {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)105081)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)60037)) {
            result[0] += -0.013439308;
          } else {
            result[0] += -0.21589331;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)1552)) {
            result[0] += -0.04809389;
          } else {
            result[0] += 0.16019432;
          }
        }
      } else {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)3889)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)15348)) {
            result[0] += 0.12485906;
          } else {
            result[0] += -0.23093326;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)299677)) {
            result[0] += -0.42783746;
          } else {
            result[0] += 0.022286218;
          }
        }
      }
    } else {
      result[0] += -0.6243624;
    }
  } else {
    if ( (data[14].missing != -1) && (data[14].fvalue < (float)31265)) {
      if ( (data[9].missing != -1) && (data[9].fvalue < (float)259658)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)889)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)302106)) {
            result[0] += 0.051358264;
          } else {
            result[0] += -0.2884319;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)125711)) {
            result[0] += 0.04851189;
          } else {
            result[0] += 0.31915456;
          }
        }
      } else {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)397)) {
          result[0] += -0.51138395;
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)89293)) {
            result[0] += 0.11253712;
          } else {
            result[0] += -0.39294267;
          }
        }
      }
    } else {
      if ( (data[0].missing != -1) && (data[0].fvalue < (float)126006)) {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)150453)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)135217)) {
            result[0] += 0.017521596;
          } else {
            result[0] += -0.17806701;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)38074)) {
            result[0] += 0.18634982;
          } else {
            result[0] += -0.049209695;
          }
        }
      } else {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)666602)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)1082)) {
            result[0] += -0.34801194;
          } else {
            result[0] += -0.086609475;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)212)) {
            result[0] += 1.2766173;
          } else {
            result[0] += -0.0143791875;
          }
        }
      }
    }
  }
  if ( (data[14].missing != -1) && (data[14].fvalue < (float)12171)) {
    if ( (data[14].missing != -1) && (data[14].fvalue < (float)5288)) {
      if ( (data[19].missing != -1) && (data[19].fvalue < (float)94033)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)86778)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)84713)) {
            result[0] += -0.03701238;
          } else {
            result[0] += 0.1352365;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)1196)) {
            result[0] += 0.18979557;
          } else {
            result[0] += -0.038826812;
          }
        }
      } else {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)5498)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)92205)) {
            result[0] += -0.09223797;
          } else {
            result[0] += 0.22591145;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)339004)) {
            result[0] += -0.3847904;
          } else {
            result[0] += 0.20550306;
          }
        }
      }
    } else {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)292)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)6678)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)360)) {
            result[0] += -0.41381183;
          } else {
            result[0] += 0.4599362;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)252839)) {
            result[0] += -0.3714623;
          } else {
            result[0] += 0.17585878;
          }
        }
      } else {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)126137)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)9309)) {
            result[0] += 0.4094616;
          } else {
            result[0] += 0.12208279;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)361615)) {
            result[0] += -0.22928552;
          } else {
            result[0] += 0.40176317;
          }
        }
      }
    }
  } else {
    if ( (data[9].missing != -1) && (data[9].fvalue < (float)890202)) {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)463148)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)179617)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)181409)) {
            result[0] += -0.025566116;
          } else {
            result[0] += -0.22124009;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)21049)) {
            result[0] += -0.054081302;
          } else {
            result[0] += -0.33929265;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)6118)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)96438)) {
            result[0] += 0.07508983;
          } else {
            result[0] += -0.42138368;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)2289)) {
            result[0] += -0.6139229;
          } else {
            result[0] += 0.25229102;
          }
        }
      }
    } else {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)1217)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)409446)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)147138)) {
            result[0] += -0.88187695;
          } else {
            result[0] += 0.15241058;
          }
        } else {
          result[0] += 0.36114076;
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)2893)) {
          result[0] += 0.0825803;
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)81234)) {
            result[0] += 0.6613543;
          } else {
            result[0] += 0.2756419;
          }
        }
      }
    }
  }
  if ( (data[10].missing != -1) && (data[10].fvalue < (float)85736)) {
    if ( (data[19].missing != -1) && (data[19].fvalue < (float)1159)) {
      if ( (data[18].missing != -1) && (data[18].fvalue < (float)337133)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)708)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)62)) {
            result[0] += -0.19147274;
          } else {
            result[0] += -0.021401403;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)365552)) {
            result[0] += 0.055365633;
          } else {
            result[0] += -0.41935048;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)72161)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)75684)) {
            result[0] += -0.9160018;
          } else {
            result[0] += -0.31567743;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)62443)) {
            result[0] += 0.1034173;
          } else {
            result[0] += -0.6791416;
          }
        }
      }
    } else {
      if ( (data[6].missing != -1) && (data[6].fvalue < (float)81901)) {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)214)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)3889)) {
            result[0] += 0.31936938;
          } else {
            result[0] += -0.015861576;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)121166)) {
            result[0] += 0.101651035;
          } else {
            result[0] += 0.27173534;
          }
        }
      } else {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)99452)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)3342)) {
            result[0] += -0.012473145;
          } else {
            result[0] += -0.28109202;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)65)) {
            result[0] += -0.23606975;
          } else {
            result[0] += 0.16689603;
          }
        }
      }
    }
  } else {
    if ( (data[22].missing != -1) && (data[22].fvalue < (float)271372)) {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)111809)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)19931)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)60490)) {
            result[0] += -0.076858856;
          } else {
            result[0] += 0.1252884;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)527264)) {
            result[0] += -0.19374685;
          } else {
            result[0] += 0.2881886;
          }
        }
      } else {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)105606)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)229774)) {
            result[0] += -0.4018186;
          } else {
            result[0] += 0.10594987;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)93210)) {
            result[0] += -0.025835147;
          } else {
            result[0] += 0.22748351;
          }
        }
      }
    } else {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)810992)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)644622)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)1108410)) {
            result[0] += -0.49452302;
          } else {
            result[0] += 0.29644373;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)5212)) {
            result[0] += 0.18398161;
          } else {
            result[0] += -0.6270633;
          }
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)128585)) {
          result[0] += -0.18765126;
        } else {
          result[0] += 0.38193485;
        }
      }
    }
  }
  if ( (data[21].missing != -1) && (data[21].fvalue < (float)1150471)) {
    if ( (data[16].missing != -1) && (data[16].fvalue < (float)273573)) {
      if ( (data[24].missing != -1) && (data[24].fvalue < (float)15830)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)6738)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)77718)) {
            result[0] += 0.030623173;
          } else {
            result[0] += -0.08492484;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)99186)) {
            result[0] += 0.13362859;
          } else {
            result[0] += -0.08361788;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)180126)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)18890)) {
            result[0] += 0.09582605;
          } else {
            result[0] += -0.08578464;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)82877)) {
            result[0] += 0.14471553;
          } else {
            result[0] += -0.099109665;
          }
        }
      }
    } else {
      if ( (data[18].missing != -1) && (data[18].fvalue < (float)128)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)239017)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)215)) {
            result[0] += -0.4939281;
          } else {
            result[0] += -0.0030498037;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)497828)) {
            result[0] += -0.752721;
          } else {
            result[0] += 0.08521007;
          }
        }
      } else {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)386)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)32256)) {
            result[0] += -0.6592328;
          } else {
            result[0] += 0.13655913;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)133222)) {
            result[0] += 0.06558017;
          } else {
            result[0] += -0.13284533;
          }
        }
      }
    }
  } else {
    if ( (data[18].missing != -1) && (data[18].fvalue < (float)299532)) {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)992579)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)1552)) {
          result[0] += 0.12251269;
        } else {
          result[0] += 0.44686386;
        }
      } else {
        result[0] += 0.04966556;
      }
    } else {
      result[0] += -0.062427502;
    }
  }
  if ( (data[4].missing != -1) && (data[4].fvalue < (float)241172)) {
    if ( (data[4].missing != -1) && (data[4].fvalue < (float)209672)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)195944)) {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)384957)) {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)9103)) {
            result[0] += 0.02191014;
          } else {
            result[0] += -0.030380178;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)389)) {
            result[0] += -0.10641497;
          } else {
            result[0] += 0.26758245;
          }
        }
      } else {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)709386)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)255869)) {
            result[0] += -0.49885884;
          } else {
            result[0] += -0.09346511;
          }
        } else {
          result[0] += 0.20177379;
        }
      }
    } else {
      if ( (data[8].missing != -1) && (data[8].fvalue < (float)17862)) {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)10350)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)61)) {
            result[0] += 0.666351;
          } else {
            result[0] += 0.049095687;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)61)) {
            result[0] += 1.0887841;
          } else {
            result[0] += 0.22940199;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)265176)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)496813)) {
            result[0] += -0.46304893;
          } else {
            result[0] += 0.069210336;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)158874)) {
            result[0] += 0.25015166;
          } else {
            result[0] += 0.8181534;
          }
        }
      }
    }
  } else {
    if ( (data[1].missing != -1) && (data[1].fvalue < (float)215018)) {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)629733)) {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)20787)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)534)) {
            result[0] += -0.3514616;
          } else {
            result[0] += 0.048353612;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)75652)) {
            result[0] += -0.5714485;
          } else {
            result[0] += -0.1678555;
          }
        }
      } else {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)58551)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)272)) {
            result[0] += 0.04688425;
          } else {
            result[0] += 0.37913036;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)183420)) {
            result[0] += -0.4165445;
          } else {
            result[0] += 0.13179682;
          }
        }
      }
    } else {
      if ( (data[9].missing != -1) && (data[9].fvalue < (float)922)) {
        result[0] += -0.3727531;
      } else {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)65300)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)22274)) {
            result[0] += 0.38205716;
          } else {
            result[0] += 1.4054606;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)10067)) {
            result[0] += -0.41679057;
          } else {
            result[0] += 0.20976107;
          }
        }
      }
    }
  }
  if ( (data[1].missing != -1) && (data[1].fvalue < (float)74117)) {
    if ( (data[12].missing != -1) && (data[12].fvalue < (float)393)) {
      if ( (data[22].missing != -1) && (data[22].fvalue < (float)138360)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)138818)) {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)217540)) {
            result[0] += -0.16040157;
          } else {
            result[0] += 0.09326543;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)333581)) {
            result[0] += 0.1568124;
          } else {
            result[0] += -0.6173657;
          }
        }
      } else {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)61842)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)88392)) {
            result[0] += 0.33141106;
          } else {
            result[0] += -0.23005866;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)135899)) {
            result[0] += 0.019646835;
          } else {
            result[0] += -0.2964807;
          }
        }
      }
    } else {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)280802)) {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)262393)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)67)) {
            result[0] += -0.003302456;
          } else {
            result[0] += 0.071616694;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)61)) {
            result[0] += 0.9979889;
          } else {
            result[0] += 0.24114402;
          }
        }
      } else {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)437876)) {
          result[0] += -0.502527;
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)222967)) {
            result[0] += -0.18571772;
          } else {
            result[0] += 0.3930854;
          }
        }
      }
    }
  } else {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)34371)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)79930)) {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)102205)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)94028)) {
            result[0] += -0.008237905;
          } else {
            result[0] += -0.51838356;
          }
        } else {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)121945)) {
            result[0] += 0.6146908;
          } else {
            result[0] += 0.079876944;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)23795)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)343786)) {
            result[0] += -0.2469174;
          } else {
            result[0] += 0.025960768;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)274090)) {
            result[0] += -0.24800558;
          } else {
            result[0] += 0.5563648;
          }
        }
      }
    } else {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)220319)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)12957)) {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)58432)) {
            result[0] += 0.0021240807;
          } else {
            result[0] += -0.3629041;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)360145)) {
            result[0] += -0.32893458;
          } else {
            result[0] += 0.027986316;
          }
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)235480)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)8285)) {
            result[0] += 0.6177016;
          } else {
            result[0] += -0.40136823;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)1354)) {
            result[0] += -0.19350626;
          } else {
            result[0] += 0.12436409;
          }
        }
      }
    }
  }
  if ( (data[22].missing != -1) && (data[22].fvalue < (float)535979)) {
    if ( (data[0].missing != -1) && (data[0].fvalue < (float)1048)) {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)204219)) {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)255965)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)58551)) {
            result[0] += -0.039514106;
          } else {
            result[0] += 0.054074157;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)124900)) {
            result[0] += -0.5130608;
          } else {
            result[0] += -0.0016904086;
          }
        }
      } else {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)108066)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)1685)) {
            result[0] += -0.7207862;
          } else {
            result[0] += -0.34209278;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)109111)) {
            result[0] += 1.0480996;
          } else {
            result[0] += -0.070257045;
          }
        }
      }
    } else {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)61174)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)55002)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)87383)) {
            result[0] += 0.0638844;
          } else {
            result[0] += 0.23635821;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)45627)) {
            result[0] += 0.10342423;
          } else {
            result[0] += -0.092996985;
          }
        }
      } else {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)23074)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)75889)) {
            result[0] += -0.2959921;
          } else {
            result[0] += 0.06242117;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)199049)) {
            result[0] += -0.17994995;
          } else {
            result[0] += 0.057661813;
          }
        }
      }
    }
  } else {
    if ( (data[11].missing != -1) && (data[11].fvalue < (float)88426)) {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)83179)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)103398)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)343219)) {
            result[0] += 0.29432103;
          } else {
            result[0] += -0.21863373;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)684)) {
            result[0] += -0.63673013;
          } else {
            result[0] += 0.035291437;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)389)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)180268)) {
            result[0] += 0.22604911;
          } else {
            result[0] += -0.3481826;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)96006)) {
            result[0] += -0.5505813;
          } else {
            result[0] += -0.040057465;
          }
        }
      }
    } else {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)105097)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)35780)) {
          result[0] += -0.77254844;
        } else {
          result[0] += -0.1308846;
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)380)) {
          result[0] += -0.37935406;
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)48261)) {
            result[0] += 0.2552673;
          } else {
            result[0] += -0.13698803;
          }
        }
      }
    }
  }
  if ( (data[11].missing != -1) && (data[11].fvalue < (float)207934)) {
    if ( (data[0].missing != -1) && (data[0].fvalue < (float)503910)) {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)265244)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)173416)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)177737)) {
            result[0] += -0.026461968;
          } else {
            result[0] += 0.07040536;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)266067)) {
            result[0] += 0.12341469;
          } else {
            result[0] += -0.29927155;
          }
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)228271)) {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)460)) {
            result[0] += 0.054605823;
          } else {
            result[0] += -0.24691926;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)1253771)) {
            result[0] += -0.6162327;
          } else {
            result[0] += 0.30740103;
          }
        }
      }
    } else {
      if ( (data[10].missing != -1) && (data[10].fvalue < (float)69617)) {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)46057)) {
          result[0] += -0.6915871;
        } else {
          result[0] += -0.33192843;
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)6547)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)14500)) {
            result[0] += 0.37972757;
          } else {
            result[0] += -0.27038953;
          }
        } else {
          result[0] += -0.4062963;
        }
      }
    }
  } else {
    if ( (data[15].missing != -1) && (data[15].fvalue < (float)1876)) {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)855)) {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)575825)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)762)) {
            result[0] += -0.5183436;
          } else {
            result[0] += -0.21577342;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)276962)) {
            result[0] += 0.7651622;
          } else {
            result[0] += -0.28267562;
          }
        }
      } else {
        if ( (data[16].missing != -1) && (data[16].fvalue < (float)288730)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)185258)) {
            result[0] += 0.2626739;
          } else {
            result[0] += -0.46489286;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)71761)) {
            result[0] += -0.5909221;
          } else {
            result[0] += 0.011317923;
          }
        }
      }
    } else {
      if ( (data[23].missing != -1) && (data[23].fvalue < (float)38074)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)211426)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)76794)) {
            result[0] += 0.28946504;
          } else {
            result[0] += 0.026024047;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)592)) {
            result[0] += -0.7285096;
          } else {
            result[0] += 0.014271273;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)124026)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)139760)) {
            result[0] += -0.31735957;
          } else {
            result[0] += 0.035695106;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)33148)) {
            result[0] += 0.26356247;
          } else {
            result[0] += -0.05388813;
          }
        }
      }
    }
  }
  if ( (data[5].missing != -1) && (data[5].fvalue < (float)65)) {
    if ( (data[15].missing != -1) && (data[15].fvalue < (float)447)) {
      if ( (data[9].missing != -1) && (data[9].fvalue < (float)41096)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)1969)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)173150)) {
            result[0] += -0.3893989;
          } else {
            result[0] += 0.1636719;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)216)) {
            result[0] += -0.1561047;
          } else {
            result[0] += 0.10647261;
          }
        }
      } else {
        if ( (data[19].missing != -1) && (data[19].fvalue < (float)1862)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)1184)) {
            result[0] += 0.41693017;
          } else {
            result[0] += -0.08991828;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)140070)) {
            result[0] += -0.29758474;
          } else {
            result[0] += 0.08308574;
          }
        }
      }
    } else {
      if ( (data[0].missing != -1) && (data[0].fvalue < (float)250202)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)174645)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)970)) {
            result[0] += -0.08221488;
          } else {
            result[0] += 0.0830609;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)1404)) {
            result[0] += 0.06239984;
          } else {
            result[0] += -0.22957423;
          }
        }
      } else {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)400469)) {
          result[0] += -0.51943654;
        } else {
          result[0] += -0.006463004;
        }
      }
    }
  } else {
    if ( (data[5].missing != -1) && (data[5].fvalue < (float)131)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)28390)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)23855)) {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)259796)) {
            result[0] += -0.420063;
          } else {
            result[0] += 0.35860997;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)9103)) {
            result[0] += 0.40420595;
          } else {
            result[0] += -0.33547565;
          }
        }
      } else {
        if ( (data[9].missing != -1) && (data[9].fvalue < (float)39671)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)62)) {
            result[0] += 1.1266086;
          } else {
            result[0] += 0.16898486;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)377)) {
            result[0] += 0.085482836;
          } else {
            result[0] += -0.3904308;
          }
        }
      }
    } else {
      if ( (data[5].missing != -1) && (data[5].fvalue < (float)224)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)252479)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)136)) {
            result[0] += 0.048724514;
          } else {
            result[0] += -0.33171675;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)61)) {
            result[0] += 0.3980194;
          } else {
            result[0] += -0.5940667;
          }
        }
      } else {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)22348)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)190)) {
            result[0] += 0.012900685;
          } else {
            result[0] += 0.10163702;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)476425)) {
            result[0] += -0.061386984;
          } else {
            result[0] += 0.2052208;
          }
        }
      }
    }
  }
  if ( (data[23].missing != -1) && (data[23].fvalue < (float)909421)) {
    if ( (data[3].missing != -1) && (data[3].fvalue < (float)1306)) {
      if ( (data[1].missing != -1) && (data[1].fvalue < (float)244411)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)201562)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)1915)) {
            result[0] += -0.049652487;
          } else {
            result[0] += 0.031230329;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)25237)) {
            result[0] += 0.17192595;
          } else {
            result[0] += -0.5348125;
          }
        }
      } else {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)72236)) {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)2307)) {
            result[0] += -0.5749784;
          } else {
            result[0] += -0.09570908;
          }
        } else {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)61620)) {
            result[0] += -0.47990805;
          } else {
            result[0] += 0.18174562;
          }
        }
      }
    } else {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)18607)) {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)1682)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)131956)) {
            result[0] += 0.11966529;
          } else {
            result[0] += 0.37940347;
          }
        } else {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)14995)) {
            result[0] += -0.16029891;
          } else {
            result[0] += 0.4682884;
          }
        }
      } else {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)31308)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)468)) {
            result[0] += -0.013102397;
          } else {
            result[0] += 0.14176564;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)389)) {
            result[0] += -0.1882702;
          } else {
            result[0] += -0.0008321063;
          }
        }
      }
    }
  } else {
    if ( (data[24].missing != -1) && (data[24].fvalue < (float)37393)) {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)303)) {
        if ( (data[0].missing != -1) && (data[0].fvalue < (float)61)) {
          result[0] += -0.6674227;
        } else {
          result[0] += 0.060863126;
        }
      } else {
        if ( (data[1].missing != -1) && (data[1].fvalue < (float)207)) {
          result[0] += 0.32304278;
        } else {
          result[0] += -0.11303438;
        }
      }
    } else {
      if ( (data[5].missing != -1) && (data[5].fvalue < (float)128902)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)992579)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)558130)) {
            result[0] += 0.37783664;
          } else {
            result[0] += -0.043974042;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)3530)) {
            result[0] += -0.23353468;
          } else {
            result[0] += 0.21194248;
          }
        }
      } else {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)78414)) {
          result[0] += -0.24848449;
        } else {
          result[0] += 0.22418244;
        }
      }
    }
  }
  if ( (data[5].missing != -1) && (data[5].fvalue < (float)107646)) {
    if ( (data[1].missing != -1) && (data[1].fvalue < (float)478932)) {
      if ( (data[5].missing != -1) && (data[5].fvalue < (float)80136)) {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)102337)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)8398)) {
            result[0] += 0.015615679;
          } else {
            result[0] += -0.090134144;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)279128)) {
            result[0] += 0.10095207;
          } else {
            result[0] += -0.45749316;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)69570)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)140011)) {
            result[0] += 0.2961992;
          } else {
            result[0] += -0.34481555;
          }
        } else {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)170738)) {
            result[0] += -0.18204941;
          } else {
            result[0] += 0.2598947;
          }
        }
      }
    } else {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)479164)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)200078)) {
          if ( (data[0].missing != -1) && (data[0].fvalue < (float)604118)) {
            result[0] += -0.6224331;
          } else {
            result[0] += -0.041272122;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)202918)) {
            result[0] += 0.3143303;
          } else {
            result[0] += -0.30297774;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)121402)) {
          result[0] += 0.055083204;
        } else {
          result[0] += 0.23806944;
        }
      }
    }
  } else {
    if ( (data[15].missing != -1) && (data[15].fvalue < (float)172251)) {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)934433)) {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)17146)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)54067)) {
            result[0] += -0.111716345;
          } else {
            result[0] += 0.16500987;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)146069)) {
            result[0] += -0.29828867;
          } else {
            result[0] += -0.07335322;
          }
        }
      } else {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)61)) {
          result[0] += -0.059009697;
        } else {
          result[0] += 0.4258038;
        }
      }
    } else {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)245283)) {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)8084)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)103225)) {
            result[0] += 0.513939;
          } else {
            result[0] += -0.32860705;
          }
        } else {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)263851)) {
            result[0] += -0.32783818;
          } else {
            result[0] += 0.287667;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)497)) {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)3405)) {
            result[0] += -0.7218274;
          } else {
            result[0] += -0.27141887;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)76794)) {
            result[0] += 0.25739345;
          } else {
            result[0] += -0.15970926;
          }
        }
      }
    }
  }
  if ( (data[20].missing != -1) && (data[20].fvalue < (float)1282683)) {
    if ( (data[24].missing != -1) && (data[24].fvalue < (float)61)) {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)266067)) {
        if ( (data[3].missing != -1) && (data[3].fvalue < (float)1172)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)384)) {
            result[0] += -0.1421551;
          } else {
            result[0] += -0.0005981191;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)268)) {
            result[0] += 0.107906245;
          } else {
            result[0] += -0.023841945;
          }
        }
      } else {
        if ( (data[12].missing != -1) && (data[12].fvalue < (float)3703)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)107646)) {
            result[0] += -0.741208;
          } else {
            result[0] += -0.23988378;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)90706)) {
            result[0] += -0.059496325;
          } else {
            result[0] += -0.40607998;
          }
        }
      }
    } else {
      if ( (data[15].missing != -1) && (data[15].fvalue < (float)218985)) {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)30542)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)88945)) {
            result[0] += 0.036897775;
          } else {
            result[0] += 0.17462674;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)252839)) {
            result[0] += -0.030089257;
          } else {
            result[0] += 0.10385827;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)148186)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)142586)) {
            result[0] += -0.31037748;
          } else {
            result[0] += -0.046265915;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)589332)) {
            result[0] += -0.034507286;
          } else {
            result[0] += 0.3944455;
          }
        }
      }
    }
  } else {
    if ( (data[5].missing != -1) && (data[5].fvalue < (float)100859)) {
      if ( (data[3].missing != -1) && (data[3].fvalue < (float)1981)) {
        result[0] += 0.4552168;
      } else {
        result[0] += 0.15239064;
      }
    } else {
      result[0] += 0.040771198;
    }
  }
  if ( (data[19].missing != -1) && (data[19].fvalue < (float)6726)) {
    if ( (data[24].missing != -1) && (data[24].fvalue < (float)200353)) {
      if ( (data[19].missing != -1) && (data[19].fvalue < (float)1862)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)355420)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)149024)) {
            result[0] += 0.007086436;
          } else {
            result[0] += 0.0922753;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)599376)) {
            result[0] += -0.55081934;
          } else {
            result[0] += 0.05571421;
          }
        }
      } else {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)129)) {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)580)) {
            result[0] += 0.11676856;
          } else {
            result[0] += -0.49907824;
          }
        } else {
          if ( (data[3].missing != -1) && (data[3].fvalue < (float)57217)) {
            result[0] += 0.290057;
          } else {
            result[0] += -0.010343186;
          }
        }
      }
    } else {
      if ( (data[21].missing != -1) && (data[21].fvalue < (float)197581)) {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)632)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)15127)) {
            result[0] += 0.15427952;
          } else {
            result[0] += -0.32395574;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)126344)) {
            result[0] += -0.50838596;
          } else {
            result[0] += 0.012661119;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)252786)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)358766)) {
            result[0] += -0.7448249;
          } else {
            result[0] += -0.18689653;
          }
        } else {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)1042)) {
            result[0] += -0.1339077;
          } else {
            result[0] += 0.28079733;
          }
        }
      }
    }
  } else {
    if ( (data[13].missing != -1) && (data[13].fvalue < (float)1292)) {
      if ( (data[19].missing != -1) && (data[19].fvalue < (float)199267)) {
        if ( (data[7].missing != -1) && (data[7].fvalue < (float)281)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)361615)) {
            result[0] += -0.0051306607;
          } else {
            result[0] += 1.1578267;
          }
        } else {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)185015)) {
            result[0] += -0.22269161;
          } else {
            result[0] += 0.07995723;
          }
        }
      } else {
        if ( (data[22].missing != -1) && (data[22].fvalue < (float)26718)) {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)7953)) {
            result[0] += -0.40087995;
          } else {
            result[0] += 0.10315466;
          }
        } else {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)793946)) {
            result[0] += -0.5114812;
          } else {
            result[0] += 0.056764092;
          }
        }
      }
    } else {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)14554)) {
        if ( (data[15].missing != -1) && (data[15].fvalue < (float)30220)) {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)225649)) {
            result[0] += 0.26712847;
          } else {
            result[0] += -0.5737403;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)4555)) {
            result[0] += 0.31956196;
          } else {
            result[0] += -0.048250135;
          }
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)456649)) {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)461)) {
            result[0] += -0.15603432;
          } else {
            result[0] += -0.014565676;
          }
        } else {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)198)) {
            result[0] += -0.19166829;
          } else {
            result[0] += 0.20759241;
          }
        }
      }
    }
  }
  if ( (data[10].missing != -1) && (data[10].fvalue < (float)18447)) {
    if ( (data[10].missing != -1) && (data[10].fvalue < (float)2080)) {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)346827)) {
        if ( (data[5].missing != -1) && (data[5].fvalue < (float)348236)) {
          if ( (data[5].missing != -1) && (data[5].fvalue < (float)321753)) {
            result[0] += 0.0055670124;
          } else {
            result[0] += 0.8628982;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)130302)) {
            result[0] += -0.52409524;
          } else {
            result[0] += 0.2132145;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)14911)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)220526)) {
            result[0] += -0.6427616;
          } else {
            result[0] += -0.15350561;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)350)) {
            result[0] += -0.6313532;
          } else {
            result[0] += -0.051136073;
          }
        }
      }
    } else {
      if ( (data[7].missing != -1) && (data[7].fvalue < (float)3791)) {
        if ( (data[2].missing != -1) && (data[2].fvalue < (float)308)) {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)119758)) {
            result[0] += 0.048339337;
          } else {
            result[0] += 0.29770347;
          }
        } else {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)323762)) {
            result[0] += -0.1892249;
          } else {
            result[0] += 0.29028258;
          }
        }
      } else {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)127)) {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)12943)) {
            result[0] += 0.12256735;
          } else {
            result[0] += -0.31541702;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)144331)) {
            result[0] += 0.33768708;
          } else {
            result[0] += 0.084161095;
          }
        }
      }
    }
  } else {
    if ( (data[21].missing != -1) && (data[21].fvalue < (float)90290)) {
      if ( (data[20].missing != -1) && (data[20].fvalue < (float)115173)) {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)62)) {
          if ( (data[12].missing != -1) && (data[12].fvalue < (float)65590)) {
            result[0] += -0.0040376624;
          } else {
            result[0] += -0.285639;
          }
        } else {
          if ( (data[16].missing != -1) && (data[16].fvalue < (float)44222)) {
            result[0] += 0.096247114;
          } else {
            result[0] += -0.11106797;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)252479)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)57512)) {
            result[0] += 0.30011663;
          } else {
            result[0] += -0.017080257;
          }
        } else {
          if ( (data[24].missing != -1) && (data[24].fvalue < (float)21546)) {
            result[0] += -0.050172087;
          } else {
            result[0] += -0.50224525;
          }
        }
      }
    } else {
      if ( (data[12].missing != -1) && (data[12].fvalue < (float)142916)) {
        if ( (data[6].missing != -1) && (data[6].fvalue < (float)174876)) {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)632)) {
            result[0] += -0.15024965;
          } else {
            result[0] += -0.38945648;
          }
        } else {
          if ( (data[7].missing != -1) && (data[7].fvalue < (float)44287)) {
            result[0] += 0.32204816;
          } else {
            result[0] += -0.13169439;
          }
        }
      } else {
        if ( (data[17].missing != -1) && (data[17].fvalue < (float)54546)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)211570)) {
            result[0] += -0.29349667;
          } else {
            result[0] += 0.24981163;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)12864)) {
            result[0] += 0.2307078;
          } else {
            result[0] += -0.01878982;
          }
        }
      }
    }
  }
  if ( (data[17].missing != -1) && (data[17].fvalue < (float)799)) {
    if ( (data[16].missing != -1) && (data[16].fvalue < (float)73677)) {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)524)) {
        if ( (data[10].missing != -1) && (data[10].fvalue < (float)322872)) {
          if ( (data[19].missing != -1) && (data[19].fvalue < (float)267282)) {
            result[0] += 0.00072239107;
          } else {
            result[0] += -0.40530726;
          }
        } else {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)178455)) {
            result[0] += -0.4908072;
          } else {
            result[0] += 0.31897745;
          }
        }
      } else {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)61)) {
          if ( (data[15].missing != -1) && (data[15].fvalue < (float)159459)) {
            result[0] += -0.37793574;
          } else {
            result[0] += -0.042292207;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)21049)) {
            result[0] += 0.058986157;
          } else {
            result[0] += -0.13868146;
          }
        }
      }
    } else {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)323290)) {
        if ( (data[11].missing != -1) && (data[11].fvalue < (float)151258)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)110247)) {
            result[0] += 0.09933209;
          } else {
            result[0] += -0.5048613;
          }
        } else {
          if ( (data[20].missing != -1) && (data[20].fvalue < (float)206272)) {
            result[0] += 0.37182632;
          } else {
            result[0] += -0.15670739;
          }
        }
      } else {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)111075)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)102181)) {
            result[0] += -0.6322176;
          } else {
            result[0] += 0.09195909;
          }
        } else {
          if ( (data[4].missing != -1) && (data[4].fvalue < (float)113484)) {
            result[0] += 0.32361087;
          } else {
            result[0] += 0.05904116;
          }
        }
      }
    }
  } else {
    if ( (data[22].missing != -1) && (data[22].fvalue < (float)234)) {
      if ( (data[16].missing != -1) && (data[16].fvalue < (float)76943)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)5694)) {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)85817)) {
            result[0] += -0.1696629;
          } else {
            result[0] += 0.17850202;
          }
        } else {
          if ( (data[8].missing != -1) && (data[8].fvalue < (float)64537)) {
            result[0] += 0.2484612;
          } else {
            result[0] += 0.006977269;
          }
        }
      } else {
        if ( (data[21].missing != -1) && (data[21].fvalue < (float)34371)) {
          if ( (data[18].missing != -1) && (data[18].fvalue < (float)337133)) {
            result[0] += -0.22496851;
          } else {
            result[0] += 0.12864883;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)256794)) {
            result[0] += 0.13204683;
          } else {
            result[0] += -0.44857463;
          }
        }
      }
    } else {
      if ( (data[24].missing != -1) && (data[24].fvalue < (float)14767)) {
        if ( (data[18].missing != -1) && (data[18].fvalue < (float)90579)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)482040)) {
            result[0] += 0.13355699;
          } else {
            result[0] += -0.22310822;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)166127)) {
            result[0] += -0.06709921;
          } else {
            result[0] += 0.21174093;
          }
        }
      } else {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)519)) {
          if ( (data[6].missing != -1) && (data[6].fvalue < (float)16166)) {
            result[0] += -0.042813692;
          } else {
            result[0] += -0.27129528;
          }
        } else {
          if ( (data[21].missing != -1) && (data[21].fvalue < (float)17345)) {
            result[0] += 0.15482087;
          } else {
            result[0] += -0.0037253245;
          }
        }
      }
    }
  }
  if ( (data[9].missing != -1) && (data[9].fvalue < (float)124)) {
    if ( (data[3].missing != -1) && (data[3].fvalue < (float)191074)) {
      if ( (data[4].missing != -1) && (data[4].fvalue < (float)67089)) {
        if ( (data[14].missing != -1) && (data[14].fvalue < (float)141680)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)332130)) {
            result[0] += 0.013023627;
          } else {
            result[0] += 0.1636096;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)180247)) {
            result[0] += -0.372462;
          } else {
            result[0] += 0.006061582;
          }
        }
      } else {
        if ( (data[8].missing != -1) && (data[8].fvalue < (float)71802)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)39954)) {
            result[0] += 0.25649175;
          } else {
            result[0] += -0.15107337;
          }
        } else {
          if ( (data[1].missing != -1) && (data[1].fvalue < (float)103786)) {
            result[0] += -0.52427447;
          } else {
            result[0] += 0.14931911;
          }
        }
      }
    } else {
      if ( (data[2].missing != -1) && (data[2].fvalue < (float)79789)) {
        if ( (data[20].missing != -1) && (data[20].fvalue < (float)257096)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)302816)) {
            result[0] += -0.70866007;
          } else {
            result[0] += -0.10019597;
          }
        } else {
          result[0] += 0.22596763;
        }
      } else {
        if ( (data[23].missing != -1) && (data[23].fvalue < (float)123849)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)144464)) {
            result[0] += -0.33166903;
          } else {
            result[0] += 0.04196429;
          }
        } else {
          if ( (data[23].missing != -1) && (data[23].fvalue < (float)214661)) {
            result[0] += 0.39550868;
          } else {
            result[0] += -0.20792098;
          }
        }
      }
    }
  } else {
    if ( (data[7].missing != -1) && (data[7].fvalue < (float)698)) {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)150547)) {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)143082)) {
          if ( (data[17].missing != -1) && (data[17].fvalue < (float)115262)) {
            result[0] += -0.11798719;
          } else {
            result[0] += 0.08345858;
          }
        } else {
          if ( (data[10].missing != -1) && (data[10].fvalue < (float)76821)) {
            result[0] += 0.062096983;
          } else {
            result[0] += 1.724033;
          }
        }
      } else {
        if ( (data[24].missing != -1) && (data[24].fvalue < (float)80308)) {
          if ( (data[2].missing != -1) && (data[2].fvalue < (float)213)) {
            result[0] += -0.5225049;
          } else {
            result[0] += -0.27775118;
          }
        } else {
          if ( (data[9].missing != -1) && (data[9].fvalue < (float)10810)) {
            result[0] += 0.2842254;
          } else {
            result[0] += -0.44487092;
          }
        }
      }
    } else {
      if ( (data[13].missing != -1) && (data[13].fvalue < (float)298)) {
        if ( (data[4].missing != -1) && (data[4].fvalue < (float)8362)) {
          if ( (data[22].missing != -1) && (data[22].fvalue < (float)293)) {
            result[0] += -0.12421276;
          } else {
            result[0] += 0.043348845;
          }
        } else {
          if ( (data[11].missing != -1) && (data[11].fvalue < (float)9768)) {
            result[0] += -0.35426423;
          } else {
            result[0] += -0.01903869;
          }
        }
      } else {
        if ( (data[13].missing != -1) && (data[13].fvalue < (float)22540)) {
          if ( (data[14].missing != -1) && (data[14].fvalue < (float)55445)) {
            result[0] += 0.025843227;
          } else {
            result[0] += 0.34996283;
          }
        } else {
          if ( (data[13].missing != -1) && (data[13].fvalue < (float)159131)) {
            result[0] += -0.08004803;
          } else {
            result[0] += 0.077611186;
          }
        }
      }
    }
  }
  
  // Apply base_scores
  result[0] += -0;
  
  // Apply postprocessor
  if (!pred_margin) { postprocess(result); }
}

void postprocess(float* result) {
  // sigmoid
  const float alpha = (float)1;
  for (size_t i = 0; i < N_TARGET * MAX_N_CLASS; ++i) {
    result[i] = (float)(1) / ((float)(1) + expf(-alpha * result[i]));
  }
}

