#ifndef KLIST_H
#define KLIST_H

// Integral cutoffs
const double CutIR = 2e-5 ;
const double CutUV = 20. ;

const double CutIRresum = 0.01 ;
const double CutUVresum = 10 ; // see M_k=10_l_lp.dat ...


/*const double klist[Nk] = { CutIR, 0.00766643,  0.0160832 ,  0.0257545 ,  0.0356934 ,0.0375, 0.04,  0.045512  ,0.048,0.051,
    0.0552959 ,  0.0652349 ,  0.0751562 ,  0.0852165 ,  0.0952608 ,
    0.105196  ,  0.115136  ,  0.125083  ,  0.135091  ,  0.145084  ,
    0.155103  ,  0.165131  ,  0.175113  ,  0.185114  ,  0.195094  ,
    0.205064  ,  0.215101  ,  0.225118  ,  0.235075  ,  0.245071  ,
    0.255076  ,  0.265072  ,  0.275063  ,  0.285057  ,  0.295074  ,
    0.305081  ,  0.315064  ,  0.325061  ,  0.335065  ,  0.345049  ,
    0.355053  ,  0.365063  ,  0.375053  ,  0.385034  ,  0.395032  ,
    0.40504   ,  0.415047  ,  0.425058  ,  0.435043  ,  0.445037  ,
    0.455044  ,  0.465038  ,  0.475048  ,  0.485041  ,  0.495024 ,
    0.52160168,  0.54817937,  0.57475705,  0.60133474,
        0.62791242,  0.65449011,  0.68106779,  0.70764547,  0.73422316,
        0.76080084,  0.78737853,  0.81395621,  0.84053389,  0.86711158,
        0.89368926,  0.92026695,  0.94684463,  0.97342232,  CutUVresum  } ; */

//const size_t Nk = 200 ;
// const double klist[Nk] = { 2e-05, 2.13633e-05, 2.28195e-05, 2.4375e-05, 2.60365e-05, 2.78112e-05, 2.9707e-05, 3.17319e-05, 3.38949e-05, 3.62053e-05, 3.86732e-05, 4.13094e-05, 4.41252e-05, 4.71329e-05, 5.03457e-05, 5.37775e-05, 5.74432e-05, 6.13588e-05, 6.55412e-05, 7.00088e-05, 7.47809e-05, 7.98783e-05, 8.53232e-05, 9.11391e-05, 9.73516e-05, 0.000103987, 0.000111076, 0.000118647, 0.000126735, 0.000135373, 0.000144601, 0.000154458, 0.000164986, 0.000176232, 0.000188245, 0.000201077, 0.000214783, 0.000229423, 0.000245062, 0.000261766, 0.000279609, 0.000298669, 0.000319027, 0.000340773, 0.000364002, 0.000388814, 0.000415317, 0.000443627, 0.000473866, 0.000506167, 0.00054067, 0.000577524, 0.000616891, 0.00065894, 0.000703857, 0.000751834, 0.000803083, 0.000857824, 0.000916297, 0.000978756, 0.00104547, 0.00111674, 0.00119286, 0.00127417, 0.00136102, 0.00145379, 0.00155289, 0.00165874, 0.00177181, 0.00189258, 0.00202159, 0.00215939, 0.00230658, 0.00246381, 0.00263175, 0.00281114, 0.00300276, 0.00320744, 0.00342608, 0.00365961, 0.00390907, 0.00417553, 0.00446015, 0.00476417, 0.00508892, 0.0054358, 0.00580633, 0.00620211, 0.00662487, 0.00707645, 0.00755881, 0.00807405, 0.00862442, 0.00921229, 0.00984024, 0.010511, 0.0112275, 0.0119928, 0.0128103, 0.0136835, 0.0146162, 0.0156125, 0.0166767, 0.0178135, 0.0190277, 0.0203247, 0.0217101, 0.02319, 0.0247707, 0.0264592, 0.0282628, 0.0301893, 0.0322471, 0.0344452, 0.0367931, 0.0393011, 0.04198, 0.0448416, 0.0478982, 0.0511631, 0.0546506, 0.0583758, 0.0623549, 0.0666053, 0.0711454, 0.075995, 0.0811751, 0.0867084, 0.0926188, 0.0989321, 0.105676, 0.112879, 0.120573, 0.128792, 0.137571, 0.146949, 0.156965, 0.167665, 0.179093, 0.191301, 0.204341, 0.21827, 0.233148, 0.24904, 0.266016, 0.284149, 0.303518, 0.324207, 0.346306, 0.369912, 0.395126, 0.42206, 0.450829, 0.48156, 0.514385, 0.549448, 0.5869, 0.626906, 0.669638, 0.715284, 0.764041, 0.816121, 0.871751, 0.931173, 0.994646, 1.06245, 1.13487, 1.21222, 1.29485, 1.38312, 1.4774, 1.5781, 1.68567, 1.80057, 1.92331, 2.05441, 2.19445, 2.34403, 2.50381, 2.67448, 2.85678, 3.05151, 3.25952, 3.4817, 3.71903, 3.97253, 4.24332, 4.53256, 4.84152, 5.17154, 5.52405, 5.90059, 6.3028, 6.73243, 7.19134, 7.68153, 8.20514, 8.76443, 9.36186, 10 } ;

const size_t Nk = 53 ;
const double klist[Nk] = { 2e-5, 5e-4, 1e-3, 4e-3, 7e-3 ,0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.43, 0.46, 0.5, 0.55, 0.62, 0.7, 0.85, 1., 1.2, 1.5, 2., 4., 10. } ;

const size_t Nkp = 1025 ; //920
const double kplist[Nkp] = { 0.01, 0.0101131, 0.0102274, 0.010343, 0.01046, 0.0105783, 0.0106979, 0.0108188, 0.0109411, 0.0110648, 0.0111899, 0.0113165, 0.0114444, 0.0115738, 0.0117047, 0.011837, 0.0119709, 0.0121062, 0.0122431, 0.0123815, 0.0125215, 0.0126631, 0.0128062, 0.012951, 0.0130975, 0.0132456, 0.0133953, 0.0135468, 0.0136999, 0.0138548, 0.0140115, 0.0141699, 0.0143301, 0.0144921, 0.014656, 0.0148217, 0.0149893, 0.0151588, 0.0153302, 0.0155035, 0.0156788, 0.0158561, 0.0160353, 0.0162166, 0.0164, 0.0165854, 0.0167729, 0.0169626, 0.0171544, 0.0173483, 0.0175445, 0.0177429, 0.0179435, 0.0181463, 0.0183515, 0.018559, 0.0187688, 0.0189811, 0.0191957, 0.0194127, 0.0196322, 0.0198542, 0.0200786, 0.0203057, 0.0205353, 0.0207674, 0.0210022, 0.0212397, 0.0214799, 0.0217227, 0.0219683, 0.0222167, 0.0224679, 0.0227219, 0.0229788, 0.0232387, 0.0235014, 0.0237671, 0.0240358, 0.0243076, 0.0245824, 0.0248604, 0.0251415, 0.0254257, 0.0257132, 0.0260039, 0.0262979, 0.0265953, 0.026896, 0.0272001, 0.0275076, 0.0278186, 0.0281332, 0.0284513, 0.0287729, 0.0290983, 0.0294273, 0.02976, 0.0300965, 0.0304368, 0.0307809, 0.0311289, 0.0314809, 0.0318368, 0.0321968, 0.0325608, 0.032929, 0.0333013, 0.0336778, 0.0340586, 0.0344437, 0.0348331, 0.0352269, 0.0356252, 0.036028, 0.0364354, 0.0368473, 0.037264, 0.0376853, 0.0381114, 0.0385423, 0.0389781, 0.0394188, 0.0398645, 0.0403152, 0.040771, 0.041232, 0.0416982, 0.0421697, 0.0426464, 0.0431286, 0.0436163, 0.0441094, 0.0446081, 0.0451125, 0.0456226, 0.0461384, 0.0466601, 0.0471876, 0.0477212, 0.0482607, 0.0488064, 0.0493582, 0.0499163, 0.0504807, 0.0510514, 0.0516286, 0.0522124, 0.0528027, 0.0533997, 0.0540035, 0.0546141, 0.0552316, 0.0558561, 0.0564876, 0.0571263, 0.0577722, 0.0584254, 0.059086, 0.059754, 0.0604296, 0.0611129, 0.0618039, 0.0625027, 0.0632093, 0.063924, 0.0646468, 0.0653777, 0.0661169, 0.0668645, 0.0676205, 0.068385, 0.0691582, 0.0699402, 0.0707309, 0.0715307, 0.0723394, 0.0731573, 0.0739845, 0.074821, 0.075667, 0.0765225, 0.0773877, 0.0782627, 0.0791476, 0.0800424, 0.0809474, 0.0818627, 0.0827883, 0.0837243, 0.0846709, 0.0856283, 0.0865964, 0.0875755, 0.0885657, 0.0895671, 0.0905798, 0.0916039, 0.0926396, 0.0936871, 0.0947464, 0.0958176, 0.096901, 0.0979966, 0.0991046, 0.100225, 0.101358, 0.102504, 0.103663, 0.104835, 0.106021, 0.107219, 0.108432, 0.109658, 0.110898, 0.112151, 0.113419, 0.114702, 0.115999, 0.11731, 0.118637, 0.119978, 0.121335, 0.122706, 0.124094, 0.125497, 0.126916, 0.128351, 0.129802, 0.13127, 0.132754, 0.134255, 0.135773, 0.137308, 0.13886, 0.14043, 0.142018, 0.143624, 0.145248, 0.14689, 0.148551, 0.15023, 0.151929, 0.153647, 0.155384, 0.157141, 0.158918, 0.160714, 0.162531, 0.164369, 0.166228, 0.168107, 0.170008, 0.17193, 0.173874, 0.17584, 0.177828, 0.179839, 0.181872, 0.183928, 0.186008, 0.188111, 0.190238, 0.192389, 0.194564, 0.196764, 0.198989, 0.201238, 0.203514, 0.205815, 0.208142, 0.210495, 0.212875, 0.215282, 0.217716, 0.220178, 0.222667, 0.225185, 0.227731, 0.230306, 0.23291, 0.235543, 0.238206, 0.2409, 0.243623, 0.246378, 0.249163, 0.251981, 0.25483, 0.257711, 0.260625, 0.263571, 0.266552, 0.269565, 0.272613, 0.275696, 0.278813, 0.281965, 0.285153, 0.288377, 0.291638, 0.294935, 0.29827, 0.301642, 0.305053, 0.308502, 0.31199, 0.315517, 0.319085, 0.322693, 0.326341, 0.330031, 0.333762, 0.337536, 0.341353, 0.345212, 0.349115, 0.353062, 0.357054, 0.361091, 0.365174, 0.369303, 0.373479, 0.377701, 0.381972, 0.386291, 0.390658, 0.395075, 0.399542, 0.404059, 0.408628, 0.413248, 0.417921, 0.422646, 0.427424, 0.432257, 0.437144, 0.442087, 0.447086, 0.452141, 0.457253, 0.462423, 0.467651, 0.472939, 0.478286, 0.483694, 0.489162, 0.494693, 0.500286, 0.505943, 0.511663, 0.517449, 0.523299, 0.529216, 0.535199, 0.541251, 0.54737, 0.553559, 0.559818, 0.566148, 0.572549, 0.579022, 0.585569, 0.59219, 0.598885, 0.605657, 0.612505, 0.61943, 0.626434, 0.633516, 0.640679, 0.647923, 0.655249, 0.662657, 0.67015, 0.677727, 0.68539, 0.693139, 0.700976, 0.708902, 0.716917, 0.725023, 0.73322, 0.74151, 0.749894, 0.758373, 0.766947, 0.775619, 0.784389, 0.793257, 0.802226, 0.811297, 0.82047, 0.829746, 0.839128, 0.848615, 0.85821, 0.867914, 0.877727, 0.887651, 0.897687, 0.907837, 0.918101, 0.928482, 0.93898, 0.949596, 0.960333, 0.971191, 0.982172, 0.993277, 1.00451, 1.01586, 1.02735, 1.03897, 1.05071, 1.06259, 1.07461, 1.08676, 1.09905, 1.11147, 1.12404, 1.13675, 1.1496, 1.1626, 1.17574, 1.18904, 1.20248, 1.21608, 1.22983, 1.24373, 1.25779, 1.27201, 1.2864, 1.30094, 1.31565, 1.33053, 1.34557, 1.36078, 1.37617, 1.39173, 1.40746, 1.42338, 1.43947, 1.45575, 1.47221, 1.48885, 1.50569, 1.52271, 1.53993, 1.55734, 1.57495, 1.59275, 1.61076, 1.62897, 1.64739, 1.66602, 1.68485, 1.7039, 1.72317, 1.74265, 1.76236, 1.78228, 1.80243, 1.82281, 1.84342, 1.86427, 1.88534, 1.90666, 1.92822, 1.95002, 1.97207, 1.99437, 2.01691, 2.03972, 2.06278, 2.0861, 2.10969, 2.13354, 2.15767, 2.18206, 2.20673, 2.23168, 2.25692, 2.28244, 2.30824, 2.33434, 2.36073, 2.38742, 2.41442, 2.44172, 2.46932, 2.49724, 2.52548, 2.55403, 2.58291, 2.61211, 2.64165, 2.67152, 2.70172, 2.73227, 2.76316, 2.7944, 2.826, 2.85795, 2.89026, 2.92294, 2.95599, 2.98941, 3.02321, 3.0574, 3.09196, 3.12692, 3.16228, 3.19803, 3.23419, 3.27076, 3.30774, 3.34514, 3.38296, 3.42121, 3.45989, 3.49901, 3.53857, 3.57858, 3.61904, 3.65996, 3.70134, 3.74319, 3.78552, 3.82832, 3.8716, 3.91538, 3.95964, 4.00441, 4.04969, 4.09548, 4.14178, 4.18861, 4.23597, 4.28387, 4.3323, 4.38129, 4.43082, 4.48092, 4.53158, 4.58282, 4.63464, 4.68704, 4.74003, 4.79363, 4.84782, 4.90264, 4.95807, 5.01413, 5.07082, 5.12815, 5.18613, 5.24477, 5.30407, 5.36404, 5.42469, 5.48603, 5.54805, 5.61078, 5.67422, 5.73838, 5.80326, 5.86887, 5.93523, 6.00234, 6.0702, 6.13883, 6.20824, 6.27844, 6.34942, 6.42121, 6.49382, 6.56724, 6.64149, 6.71658, 6.79253, 6.86932, 6.94699, 7.02554, 7.10497, 7.18531, 7.26655, 7.34871, 7.4318, 7.51582, 7.6008, 7.68674, 7.77365, 7.86154, 7.95043, 8.04032, 8.13123, 8.22317, 8.31614, 8.41017, 8.50526, 8.60142, 8.69868, 8.79703, 8.89649, 8.99708, 9.09881, 9.20168, 9.30572, 9.41094, 9.51734, 9.62495, 9.73377, 9.84383, 9.95513, 10.0677, 10.1815, 10.2966, 10.4131, 10.5308, 10.6499, 10.7703, 10.892, 11.0152, 11.1397, 11.2657, 11.3931, 11.5219, 11.6522, 11.7839, 11.9171, 12.0519, 12.1881, 12.3259, 12.4653, 12.6063, 12.7488, 12.8929, 13.0387, 13.1861, 13.3352, 13.486, 13.6385, 13.7927, 13.9486, 14.1063, 14.2658, 14.4271, 14.5902, 14.7552, 14.922, 15.0908, 15.2614, 15.4339, 15.6084, 15.7849, 15.9634, 16.1439, 16.3264, 16.511, 16.6977, 16.8865, 17.0774, 17.2705, 17.4658, 17.6632, 17.8629, 18.0649, 18.2692, 18.4757, 18.6846, 18.8959, 19.1095, 19.3256, 19.5441, 19.7651, 19.9885, 20.2145, 20.4431, 20.6742, 20.908, 21.1444, 21.3835, 21.6252, 21.8697, 22.117, 22.3671, 22.62, 22.8757, 23.1344, 23.3959, 23.6605, 23.928, 24.1985, 24.4721, 24.7488, 25.0287, 25.3116, 25.5978, 25.8873, 26.1799, 26.476, 26.7753, 27.078, 27.3842, 27.6938, 28.0069, 28.3236, 28.6438, 28.9677, 29.2952, 29.6265, 29.9614, 30.3002, 30.6428, 30.9892, 31.3396, 31.694, 32.0523, 32.4147, 32.7812, 33.1519, 33.5267, 33.9058, 34.2891, 34.6768, 35.0689, 35.4654, 35.8664, 36.2719, 36.682, 37.0968, 37.5162, 37.9404, 38.3693, 38.8032, 39.2419, 39.6856, 40.1343, 40.5881, 41.047, 41.5111, 41.9804, 42.4551, 42.9351, 43.4205, 43.9115, 44.408, 44.9101, 45.4178, 45.9314, 46.4507, 46.9759, 47.507, 48.0442, 48.5874, 49.1367, 49.6923, 50.2541, 50.8223, 51.397, 51.9781, 52.5658, 53.1601, 53.7612, 54.369, 54.9838, 55.6054, 56.2341, 56.8699, 57.5129, 58.1632, 58.8208, 59.4859, 60.1585, 60.8387, 61.5265, 62.2222, 62.9257, 63.6372, 64.3567, 65.0843, 65.8202, 66.5644, 67.317, 68.0782, 68.8479, 69.6263, 70.4136, 71.2097, 72.0148, 72.8291, 73.6525, 74.4853, 75.3274, 76.1791, 77.0404, 77.9115, 78.7924, 79.6833, 80.5842, 81.4953, 82.4168, 83.3486, 84.291, 85.244, 86.2079, 87.1826, 88.1683, 89.1652, 90.1733, 91.1929, 92.224, 93.2667, 94.3212, 95.3877, 96.4662, 97.5569, 98.6599, 99.7754, 100.904, 102.044, 103.198, 104.365, 105.545, 106.738, 107.945, 109.166, 110.4, 111.648, 112.911, 114.187, 115.478, 116.784, 118.104, 119.44, 120.79, 122.156, 123.537, 124.934, 126.346, 127.775, 129.22, 130.681, 132.158, 133.652, 135.163, 136.692, 138.237, 139.8, 141.381, 142.979, 144.596, 146.231, 147.884, 149.556, 151.247, 152.957, 154.687, 156.436, 158.204, 159.993, 161.802, 163.632, 165.482, 167.353, 169.245, 171.158, 173.094, 175.051, 177.03, 179.032, 181.056, 183.103, 185.173, 187.267, 189.384, 191.525, 193.691, 195.881, 198.096, 200.335, 202.601, 204.891, 207.208, 209.551, 211.92, 214.316, 216.739, 219.19, 221.668, 224.174, 226.709, 229.272, 231.865, 234.486, 237.137, 239.819, 242.53, 245.272, 248.045, 250.85, 253.686, 256.555, 259.455, 262.389, 265.356, 268.356, 271.39, 274.458, 277.562, 280.7, 283.874, 287.083, 290.329, 293.612, 296.931, 300.289, 303.684, 307.118, 310.59, 314.102, 317.653, 321.245, 324.877, 328.55, 332.265, 336.022, 339.821, 343.663, 347.549, 351.478, 355.452, 359.471, 363.536, 367.646, 371.803, 376.006, 380.258, 384.557, 388.905, 393.302, 397.749, 402.246, 406.794, 411.394, 416.045, 420.749, 425.507, 430.318, 435.183, 440.103, 445.079, 450.112, 455.201, 460.348, 465.553, 470.816, 476.14, 481.523, 486.968, 492.473, 498.042, 503.673, 509.368, 515.127, 520.951, 526.841, 532.798, 538.822, 544.914, 551.075, 557.306, 563.607, 569.98, 576.424, 582.942, 589.533, 596.198, 602.939, 609.756, 616.65, 623.623, 630.674, 637.804, 645.016, 652.309, 659.684, 667.143, 674.686, 682.314, 690.029, 697.831, 705.721, 713.7, 721.769, 729.93, 738.183, 746.529, 754.97, 763.506, 772.139, 780.869, 789.698, 798.627, 807.656, 816.788, 826.023, 835.363, 844.808, 854.359, 864.019, 873.788, 883.668, 893.659, 903.763, 913.982, 924.316, 934.766, 945.335, 956.024, 966.833, 977.765, 988.82, 1000 } ;

const double ps1D_errrel[Nk] = { 0.05, 0.044,0.016,0.010,0.007,0.006,0.005,0.004,0.004,0.003,0.003,0.003,0.003,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001, 0.001 } ;
const double ps1D_inverrelmax = 1000. ;

#endif