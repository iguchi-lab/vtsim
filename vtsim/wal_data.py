wall_basic = {
    '木造_床_無断熱':           [[[['合板',                                 0.120 ]], 1.00]],
    '木造_間仕切壁_2重中空':    [[[['せっこうボード',                       0.0125 ],
                                   ['中空層',                               0.100 ],
                                   ['せっこうボード',                       0.0125 ]], 1.00]],
    '木造_天井_無断熱':         [[[['せっこうボード',                       0.0125]], 1.00]],
}

wall_kameido= {
     '亀戸_外壁':   [[[['せっこうボード',                       0.0125],
                       ['中空層',                               0.040],
                       ['硬質ウレタンフォーム保温板1種2号',     0.030],
                       ['ALC',                                  0.100]],    1.00]],
    '亀戸_戸境壁':  [[[['せっこうボード',                       0.0095],
                       ['せっこうボード',                       0.021],
                       ['中空層',                               0.0125],
                       ['グラスウール断熱材24K相当',            0.050],
                       ['中空層',                               0.0125],
                       ['せっこうボード',                       0.021],
                       ['せっこうボード',                       0.0095]],   1.00]],
    '亀戸_床板':    [[[['合板',                                 0.012],
                       ['せっこうボード',                       0.0095],
                       ['パーティクルボード',                   0.020]],    1.00]],
    '亀戸_スラブ':  [[[['せっこうボード',                       0.0095],
                       ['中空層',                               0.1505],
                       ['硬質ウレタンフォーム保温板1種2号',     0.020],
                       ['コンクリート',                         0.350]],    1.00]
    ]
}

wall_FPJ = {
    'FPJ_RC_基礎壁_内断熱50+30':            [[[['押出法ポリスチレンフォーム3種',        0.030],
                                               ['コンクリート',                         0.150],
                                               ['押出法ポリスチレンフォーム3種',        0.050],
                                               ['セメント・モルタル',                   0.150]],    1.00]],
    'FPJ_RC_床_内断熱30':                   [[[['押出法ポリスチレンフォーム3種',        0.030],
                                               ['コンクリート',                         0.150 ]],   1.00]],
    'FPJ_木造（在来）_外壁_充填+外断熱':    [[[['せっこうボード',                       0.0125],
                                               ['吹込用グラスウール断熱材2種35K相当',   0.120],
                                               ['通気層',                               0.180],
                                               ['木片セメント板',                       0.015],
                                               ['硬質ウレタンフォーム保温板1種2号',     0.030]],    0.83],
                                             [[['せっこうボード',                       0.0125],
                                               ['天然木材1類(桧、杉、えぞ松等)',        0.120],
                                               ['通気層',                               0.180],
                                               ['木片セメント板',                       0.015],
                                               ['硬質ウレタンフォーム保温板1種2号',     0.030]],    0.17]],
    'FPJ_木造（在来）_外壁_充填':           [[[['せっこうボード',                       0.0125],
                                               ['吹込用グラスウール断熱材2種35K相当',   0.120],
                                               ['通気層',                               0.180],
                                               ['木片セメント板',                       0.015]],    0.83],
                                             [[['せっこうボード',                       0.0125],
                                               ['天然木材1類(桧、杉、えぞ松等)',        0.120],
                                               ['通気層',                               0.180],
                                               ['木片セメント板',                       0.015]],    0.17]],        
    'FPJ_木造_天井_吹込':                   [[[['せっこうボード',                       0.0125],
                                               ['吹込用グラスウール断熱材2種35K相当',   0.200]],    1.00]],
    'FPJ_木造_間仕切壁':                    [[[['押出法ポリスチレンフォーム3種',        0.050]],    0.83],
                                             [[['天然木材1類(桧、杉、えぞ松等)',        0.120],
                                               ['押出法ポリスチレンフォーム3種',        0.050]],    0.17]],
    'FPJ_木造_屋根_充填':                   [[[['吹込用グラスウール断熱材2種35K相当',   0.200],
                                               ['合板',                                 0.012],
                                               ['通気層',                               0.030],
                                               ['合板',                                 0.012]],    0.86],
                                             [[['天然木材1類(桧、杉、えぞ松等)',        0.200],
                                               ['合板',                                 0.012],
                                               ['通気層',                               0.030],
                                               ['合板',                                 0.012]],    0.14]],
    'FPJ_RC_床_内断熱':                     [[[['押出法ポリスチレンフォーム3種',        0.015 ]],   1.00]]
}

wall_okayama = {
    '岡山_外壁':                    [[[['せっこうボード',                   0.013],
                                       ['中空層',                           0.055],
                                       ['硬質ウレタン保温材',               0.050],
                                       ['ダイライトMS',                     0.009]],    0.83],
                                     [[['せっこうボード',                   0.013],
                                       ['中空層',                           0.055],
                                       ['天然木材1類(桧、杉、えぞ松等)',    0.050],
                                       ['ダイライトMS',                     0.009]],    0.17]],
    '岡山_床(上部)':                [[[['合板',                             0.012],
                                       ['合板',                             0.012],
                                       ['ふく射パネル_PP',                  0.020]],    1.00]],
    '岡山_床(下部)':                [[[['ふく射パネル_PP',                  0.020],
                                       ['ふく射パネル_PB',                  0.020]],    1.00]],
    '岡山_床(最下部)':              [[[['合板',                             0.020],
                                       ['押出法ポリスチレンフォーム2種',    0.030]],    0.85],
                                     [[['合板',                             0.020],
                                       ['合板',                             0.030]],    0.15]],
    '岡山_天井':                    [[[['合板',                             0.024],
                                       ['遮音マット',                       0.009],
                                       ['合板',                             0.012]],    1.00]],
    '岡山_人工気象室_外壁':         [[[['合板',                             0.012],
                                       ['グラスウール断熱材10K',            0.100],
                                       ['ダイライトMS',                     0.009]],    0.82],
                                     [[['合板',                             0.012],
                                      ['グラスウール断熱材10K',             0.100],
                                      ['天然木材1類(桧、杉、えぞ松等)',     0.105],
                                      ['ダイライトMS',                      0.009]],    0.18]],
    '岡山_人工気象室_屋根':         [[[['XPS1種A',                          0.030],
                                       ['天然木材1類(桧、杉、えぞ松等)',    0.030],
                                       ['合板',                             0.009]],    0.86],
                                     [[['XPS1種A',                          0.012],
                                       ['グラスウール断熱材10K',            0.030],
                                       ['天然木材1類(桧、杉、えぞ松等)',    0.030],
                                       ['合板',                             0.009]],    0.14]],
    '岡山_人工気象室_床':           [[[['合板',                             0.024],
                                       ['XPS3種',                           0.030]],    0.81],
                                     [[['合板',                             0.024],
                                       ['天然木材1類(桧、杉、えぞ松等)',    0.030]],    0.19]],
    '岡山_人工気象室_パネル_上':    [[[['パーティクルボード',               0.020]],    1.00]],
    '岡山_人工気象室_パネル_下':    [[[['樹脂パネル',                       0.00314]],  1.00]],
}

wall_oosaka = {
    '大阪_天井':            [[[['高性能グラスウール断熱材16K相当',      0.200 ],
                               ['せっこうボード',                       0.0095]],   1.00]],
    '大阪_壁':              [[[['天然木材1類(桧、杉、えぞ松等)',        0.105 ],
                               ['合板',                                 0.012 ]],   0.17],
                             [[['高性能グラスウール断熱材16K相当',      0.105 ],
                               ['合板',                                 0.012 ]],   0.83]],
    '大阪_床':              [[[['合板',                                 0.024 ],
                               ['合板',                                 0.012 ]],   1.00]],
    '大阪_畳床':            [[[['合板',                                 0.024 ],
                               ['畳',                                   0.060 ]],   1.00]],
    '大阪_床下外壁':        [[[['A種フェノールフォーム保温板1種2号',    0.050 ],
                               ['コンクリート',                         0.150 ]],   1.00]],  
    '大阪_床下外壁_断熱無': [[[['コンクリート',                         0.150 ]],   1.00]],  
    '大阪_床下内壁':        [[[['コンクリート',                         0.150 ]],   1.00]],  
    '大阪_床下床':          [[[['A種フェノールフォーム保温板1種2号',    0.035 ],
                               ['コンクリート',                         0.150 ]],   1.00]]  
}