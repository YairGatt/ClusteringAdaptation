#!/usr/bin/env python

"""
Metabolic plotting functions for metabolic_clustering
"""

import sys, os
#from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from pandas import DataFrame
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from itertools import combinations
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
import itertools
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def compute_jaccard_index(set_1, set_2):
	set_1 = set(set_1)
	set_2 = set(set_2)
	n = len(set_1.intersection(set_2))
	try:
		return n / float(len(set_1) + len(set_2) - n)
	except ZeroDivisionError:
		return 1

def clustering_pathways_by_jaccard(clustering_dict, cluster_outfile):
	#cluster patients before continuing
	raw_linkage = [1-compute_jaccard_index(i[0],i[1]) for i in list(combinations(clustering_dict.values(),2))]
	#raw_linkage = [1/compute_jaccard_index(i[0],i[1])for i in list(combinations(patient_dict.values(),2))]
	Z = linkage(raw_linkage, 'average')
	c,coph = cophenet(Z,raw_linkage)
	# calculate full dendrogram
	fig = plt.figure(figsize=(25, 22))
	#plt.figure()
	plt.title('Hierarchical Clustering Dendrogram', fontsize=20)
	plt.xlabel('Patient', fontsize=20)
	plt.ylabel('Distance (1-Jaccard index)', fontsize=20)
	dendrogram(Z, leaf_rotation=90., leaf_font_size=16., labels=clustering_dict.keys())
	#plt.gcf()
	#dendrogram(Z, truncate_mode='lastp', p=12, leaf_rotation=90., leaf_font_size=12., show_contracted=True)
	#cluster_outfile = "_cluster.".join("test_stuff.pdf".rsplit(".",1))
	plt.savefig(cluster_outfile)

def plot_pca_no_labels(outfile,names, array, title="Default"):
        """
        Plot PCA plot in 2D with no labels
        """
        #perform pca
        pca = PCA(n_components=2)
        pca.fit(array)
        explained = pca.explained_variance_ratio_
        #reduce dimensions
        X = pca.transform(array)
        #plot
        plt.clf()
        plt.cla()
        plt.close()
        fig = plt.figure(1, figsize=(25, 22))
        #scatter plot
        plt.scatter(X[:,0],X[:,1])
        #label axes
        plt.xlabel('First principal component: %.3f' % explained[0], fontsize=20)
        plt.ylabel('Second principal component: %.3f' % explained[1], fontsize=20)
        #add title
        plt.title('2D Principal Component Analysis ' + title, fontsize=20)
        #add text
        for n, name in enumerate(names): plt.text(X[n,0], X[n,1], name, fontsize=10)
        #save or show
        plt.savefig(outfile)
        return names, X, explained


def plot_pca(outfile, names, array, labels, probabilities, components=2, title="Default"):
	"""
	Plot PCA plot in 2D or 3D
	"""
	print("Num_clusters: %s" % (max(labels)+1))
	#color map
	#LABEL_COLOR_MAP = {0:"b",1:"r",2:"g",3:"yellow",4:"c",5:"m",6:"k",7:"gray",8:"gold",9:"magenta",10:"y",11:"teal",12:"orange",13:"maroon",14:"skyblue",15:"pink",16:"navy",17:"seagreen",18:"beige",19:"salmon",20:"indigo"}
	LABEL_COLOR_MAP = {0: '#FFFF00', 1: '#1CE6FF', 2: '#FF34FF', 3: '#FF4A46', 4: '#008941', 5: '#006FA6', 6: '#A30059', 7: '#FFDBE5', 8: '#7A4900', 9: '#0000A6', 10: '#63FFAC', 11: '#B79762', 12: '#004D43', 13: '#8FB0FF', 14: '#997D87', 15: '#5A0007', 16: '#809693', 17: '#FEFFE6', 18: '#1B4400', 19: '#4FC601', 20: '#3B5DFF', 21: '#4A3B53', 22: '#FF2F80', 23: '#61615A', 24: '#BA0900', 25: '#6B7900', 26: '#00C2A0', 27: '#FFAA92', 28: '#FF90C9', 29: '#B903AA', 30: '#D16100', 31: '#DDEFFF', 32: '#000035', 33: '#7B4F4B', 34: '#A1C299', 35: '#300018', 36: '#0AA6D8', 37: '#013349', 38: '#00846F', 39: '#372101', 40: '#FFB500', 41: '#C2FFED', 42: '#A079BF', 43: '#CC0744', 44: '#C0B9B2', 45: '#C2FF99', 46: '#001E09', 47: '#00489C', 48: '#6F0062', 49: '#0CBD66', 50: '#EEC3FF', 51: '#456D75', 52: '#B77B68', 53: '#7A87A1', 54: '#788D66', 55: '#885578', 56: '#FAD09F', 57: '#FF8A9A', 58: '#D157A0', 59: '#BEC459', 60: '#456648', 61: '#0086ED', 62: '#886F4C', 63: '#34362D', 64: '#B4A8BD', 65: '#00A6AA', 66: '#452C2C', 67: '#636375', 68: '#A3C8C9', 69: '#FF913F', 70: '#938A81', 71: '#575329', 72: '#00FECF', 73: '#B05B6F', 74: '#8CD0FF', 75: '#3B9700', 76: '#04F757', 77: '#C8A1A1', 78: '#1E6E00', 79: '#7900D7', 80: '#A77500', 81: '#6367A9', 82: '#A05837', 83: '#6B002C', 84: '#772600', 85: '#D790FF', 86: '#9B9700', 87: '#549E79', 88: '#FFF69F', 89: '#201625', 90: '#72418F', 91: '#BC23FF', 92: '#99ADC0', 93: '#3A2465', 94: '#922329', 95: '#5B4534', 96: '#FDE8DC', 97: '#404E55', 98: '#0089A3', 99: '#CB7E98', 100: '#A4E804', 101: '#324E72', 102: '#6A3A4C', 103: '#83AB58', 104: '#001C1E', 105: '#D1F7CE', 106: '#004B28', 107: '#C8D0F6', 108: '#A3A489', 109: '#806C66', 110: '#222800', 111: '#BF5650', 112: '#E83000', 113: '#66796D', 114: '#DA007C', 115: '#FF1A59', 116: '#8ADBB4', 117: '#1E0200', 118: '#5B4E51', 119: '#C895C5', 120: '#320033', 121: '#FF6832', 122: '#66E1D3', 123: '#CFCDAC', 124: '#D0AC94', 125: '#7ED379', 126: '#012C58', 127: '#7A7BFF', 128: '#D68E01', 129: '#353339', 130: '#78AFA1', 131: '#FEB2C6', 132: '#75797C', 133: '#837393', 134: '#943A4D', 135: '#B5F4FF', 136: '#D2DCD5', 137: '#9556BD', 138: '#6A714A', 139: '#001325', 140: '#02525F', 141: '#0AA3F7', 142: '#E98176', 143: '#DBD5DD', 144: '#5EBCD1', 145: '#3D4F44', 146: '#7E6405', 147: '#02684E', 148: '#962B75', 149: '#8D8546', 150: '#9695C5', 151: '#E773CE', 152: '#D86A78', 153: '#3E89BE', 154: '#CA834E', 155: '#518A87', 156: '#5B113C', 157: '#55813B', 158: '#E704C4', 159: '#00005F', 160: '#A97399', 161: '#4B8160', 162: '#59738A', 163: '#FF5DA7', 164: '#F7C9BF', 165: '#643127', 166: '#513A01', 167: '#6B94AA', 168: '#51A058', 169: '#A45B02', 170: '#1D1702', 171: '#E20027', 172: '#E7AB63', 173: '#4C6001', 174: '#9C6966', 175: '#64547B', 176: '#97979E', 177: '#006A66', 178: '#391406', 179: '#F4D749', 180: '#0045D2', 181: '#006C31', 182: '#DDB6D0', 183: '#7C6571', 184: '#9FB2A4', 185: '#00D891', 186: '#15A08A', 187: '#BC65E9', 188: '#FFFFFE', 189: '#C6DC99', 190: '#203B3C', 191: '#671190', 192: '#6B3A64', 193: '#F5E1FF', 194: '#FFA0F2', 195: '#CCAA35', 196: '#374527', 197: '#8BB400', 198: '#797868', 199: '#C6005A', 200: '#3B000A', 201: '#C86240', 202: '#29607C', 203: '#402334', 204: '#7D5A44', 205: '#CCB87C', 206: '#B88183', 207: '#AA5199', 208: '#B5D6C3', 209: '#A38469', 210: '#9F94F0', 211: '#A74571', 212: '#B894A6', 213: '#71BB8C', 214: '#00B433', 215: '#789EC9', 216: '#6D80BA', 217: '#953F00', 218: '#5EFF03', 219: '#E4FFFC', 220: '#1BE177', 221: '#BCB1E5', 222: '#76912F', 223: '#003109', 224: '#0060CD', 225: '#D20096', 226: '#895563', 227: '#29201D', 228: '#5B3213', 229: '#A76F42', 230: '#89412E', 231: '#1A3A2A', 232: '#494B5A', 233: '#A88C85', 234: '#F4ABAA', 235: '#A3F3AB', 236: '#00C6C8', 237: '#EA8B66', 238: '#958A9F', 239: '#BDC9D2', 240: '#9FA064', 241: '#BE4700', 242: '#658188', 243: '#83A485', 244: '#453C23', 245: '#47675D', 246: '#3A3F00', 247: '#061203', 248: '#DFFB71', 249: '#868E7E', 250: '#98D058', 251: '#6C8F7D', 252: '#D7BFC2', 253: '#3C3E6E', 254: '#D83D66', 255: '#2F5D9B', 256: '#6C5E46', 257: '#D25B88', 258: '#5B656C', 259: '#00B57F', 260: '#545C46', 261: '#866097', 262: '#365D25', 263: '#252F99', 264: '#00CCFF', 265: '#674E60', 266: '#FC009C', 267: '#92896B', 268: '#1E2324', 269: '#DEC9B2', 270: '#9D4948', 271: '#85ABB4', 272: '#342142', 273: '#D09685', 274: '#A4ACAC', 275: '#00FFFF', 276: '#AE9C86', 277: '#742A33', 278: '#0E72C5', 279: '#AFD8EC', 280: '#C064B9', 281: '#91028C', 282: '#FEEDBF', 283: '#FFB789', 284: '#9CB8E4', 285: '#AFFFD1', 286: '#2A364C', 287: '#4F4A43', 288: '#647095', 289: '#34BBFF', 290: '#807781', 291: '#920003', 292: '#B3A5A7', 293: '#018615', 294: '#F1FFC8', 295: '#976F5C', 296: '#FF3BC1', 297: '#FF5F6B', 298: '#077D84', 299: '#F56D93', 300: '#5771DA', 301: '#4E1E2A', 302: '#830055', 303: '#02D346', 304: '#BE452D', 305: '#00905E', 306: '#BE0028', 307: '#6E96E3', 308: '#007699', 309: '#FEC96D', 310: '#9C6A7D', 311: '#3FA1B8', 312: '#893DE3', 313: '#79B4D6', 314: '#7FD4D9', 315: '#6751BB', 316: '#B28D2D', 317: '#E27A05', 318: '#DD9CB8', 319: '#AABC7A', 320: '#980034', 321: '#561A02', 322: '#8F7F00', 323: '#635000', 324: '#CD7DAE', 325: '#8A5E2D', 326: '#FFB3E1', 327: '#6B6466', 328: '#C6D300', 329: '#0100E2', 330: '#88EC69', 331: '#8FCCBE', 332: '#21001C', 333: '#511F4D', 334: '#E3F6E3', 335: '#FF8EB1', 336: '#6B4F29', 337: '#A37F46', 338: '#6A5950', 339: '#1F2A1A', 340: '#04784D', 341: '#101835', 342: '#E6E0D0', 343: '#FF74FE', 344: '#00A45F', 345: '#8F5DF8', 346: '#4B0059', 347: '#412F23', 348: '#D8939E', 349: '#DB9D72', 350: '#604143', 351: '#B5BACE', 352: '#989EB7', 353: '#D2C4DB', 354: '#A587AF', 355: '#77D796', 356: '#7F8C94', 357: '#FF9B03', 358: '#555196', 359: '#31DDAE', 360: '#74B671', 361: '#802647', 362: '#2A373F', 363: '#014A68', 364: '#696628', 365: '#4C7B6D', 366: '#002C27', 367: '#7A4522', 368: '#3B5859', 369: '#E5D381', 370: '#FFF3FF', 371: '#679FA0', 372: '#261300', 373: '#2C5742', 374: '#9131AF', 375: '#AF5D88', 376: '#C7706A', 377: '#61AB1F', 378: '#8CF2D4', 379: '#C5D9B8', 380: '#9FFFFB', 381: '#BF45CC', 382: '#493941', 383: '#863B60', 384: '#B90076', 385: '#003177', 386: '#C582D2', 387: '#C1B394', 388: '#602B70', 389: '#887868', 390: '#BABFB0', 391: '#030012', 392: '#D1ACFE', 393: '#7FDEFE', 394: '#4B5C71', 395: '#A3A097', 396: '#E66D53', 397: '#637B5D', 398: '#92BEA5', 399: '#00F8B3', 400: '#BEDDFF', 401: '#3DB5A7', 402: '#DD3248', 403: '#B6E4DE', 404: '#427745', 405: '#598C5A', 406: '#B94C59', 407: '#8181D5', 408: '#94888B', 409: '#FED6BD', 410: '#536D31', 411: '#6EFF92', 412: '#E4E8FF', 413: '#20E200', 414: '#FFD0F2', 415: '#4C83A1', 416: '#BD7322', 417: '#915C4E', 418: '#8C4787', 419: '#025117', 420: '#A2AA45', 421: '#2D1B21', 422: '#A9DDB0', 423: '#FF4F78', 424: '#528500', 425: '#009A2E', 426: '#17FCE4', 427: '#71555A', 428: '#525D82', 429: '#00195A', 430: '#967874', 431: '#555558', 432: '#0B212C', 433: '#1E202B', 434: '#EFBFC4', 435: '#6F9755', 436: '#6F7586', 437: '#501D1D', 438: '#372D00', 439: '#741D16', 440: '#5EB393', 441: '#B5B400', 442: '#DD4A38', 443: '#363DFF', 444: '#AD6552', 445: '#6635AF', 446: '#836BBA', 447: '#98AA7F', 448: '#464836', 449: '#322C3E', 450: '#7CB9BA', 451: '#5B6965', 452: '#707D3D', 453: '#7A001D', 454: '#6E4636', 455: '#443A38', 456: '#AE81FF', 457: '#489079', 458: '#897334', 459: '#009087', 460: '#DA713C', 461: '#361618', 462: '#FF6F01', 463: '#006679', 464: '#370E77', 465: '#4B3A83', 466: '#C9E2E6', 467: '#C44170', 468: '#FF4526', 469: '#73BE54', 470: '#C4DF72', 471: '#ADFF60', 472: '#00447D', 473: '#DCCEC9', 474: '#BD9479', 475: '#656E5B', 476: '#EC5200', 477: '#FF6EC2', 478: '#7A617E', 479: '#DDAEA2', 480: '#77837F', 481: '#A53327', 482: '#608EFF', 483: '#B599D7', 484: '#A50149', 485: '#4E0025', 486: '#C9B1A9', 487: '#03919A', 488: '#1B2A25', 489: '#E500F1', 490: '#982E0B', 491: '#B67180', 492: '#E05859', 493: '#006039', 494: '#578F9B', 495: '#305230', 496: '#CE934C', 497: '#B3C2BE', 498: '#C0BAC0', 499: '#B506D3', 500: '#170C10', 501: '#4C534F', 502: '#224451', 503: '#3E4141', 504: '#78726D', 505: '#B6602B', 506: '#200441', 507: '#DDB588', 508: '#497200', 509: '#C5AAB6', 510: '#033C61', 511: '#71B2F5', 512: '#A9E088', 513: '#4979B0', 514: '#A2C3DF', 515: '#784149', 516: '#2D2B17', 517: '#3E0E2F', 518: '#57344C', 519: '#0091BE', 520: '#E451D1', 521: '#4B4B6A', 522: '#5C011A', 523: '#7C8060', 524: '#FF9491', 525: '#4C325D', 526: '#005C8B', 527: '#E5FDA4', 528: '#68D1B6', 529: '#032641', 530: '#140023', 531: '#8683A9', 532: '#CFFF00', 533: '#A72C3E', 534: '#34475A', 535: '#B1BB9A', 536: '#B4A04F', 537: '#8D918E', 538: '#A168A6', 539: '#813D3A', 540: '#425218', 541: '#DA8386', 542: '#776133', 543: '#563930', 544: '#8498AE', 545: '#90C1D3', 546: '#B5666B', 547: '#9B585E', 548: '#856465', 549: '#AD7C90', 550: '#E2BC00', 551: '#E3AAE0', 552: '#B2C2FE', 553: '#FD0039', 554: '#009B75', 555: '#FFF46D', 556: '#E87EAC', 557: '#DFE3E6', 558: '#848590', 559: '#AA9297', 560: '#83A193', 561: '#577977', 562: '#3E7158', 563: '#C64289', 564: '#EA0072', 565: '#C4A8CB', 566: '#55C899', 567: '#E78FCF', 568: '#004547', 569: '#F6E2E3', 570: '#966716', 571: '#378FDB', 572: '#435E6A', 573: '#DA0004', 574: '#1B000F', 575: '#5B9C8F', 576: '#6E2B52', 577: '#011115', 578: '#E3E8C4', 579: '#AE3B85', 580: '#EA1CA9', 581: '#FF9E6B', 582: '#457D8B', 583: '#92678B', 584: '#00CDBB', 585: '#9CCC04', 586: '#002E38', 587: '#96C57F', 588: '#CFF6B4', 589: '#492818', 590: '#766E52', 591: '#20370E', 592: '#E3D19F', 593: '#2E3C30', 594: '#B2EACE', 595: '#F3BDA4', 596: '#A24E3D', 597: '#976FD9', 598: '#8C9FA8', 599: '#7C2B73', 600: '#4E5F37', 601: '#5D5462', 602: '#90956F', 603: '#6AA776', 604: '#DBCBF6', 605: '#DA71FF', 606: '#987C95', 607: '#52323C', 608: '#BB3C42', 609: '#584D39', 610: '#4FC15F', 611: '#A2B9C1', 612: '#79DB21', 613: '#1D5958', 614: '#BD744E', 615: '#160B00', 616: '#20221A', 617: '#6B8295', 618: '#00E0E4', 619: '#102401', 620: '#1B782A', 621: '#DAA9B5', 622: '#B0415D', 623: '#859253', 624: '#97A094', 625: '#06E3C4', 626: '#47688C', 627: '#7C6755', 628: '#075C00', 629: '#7560D5', 630: '#7D9F00', 631: '#C36D96', 632: '#4D913E', 633: '#5F4276', 634: '#FCE4C8', 635: '#303052', 636: '#4F381B', 637: '#E5A532', 638: '#706690', 639: '#AA9A92', 640: '#237363', 641: '#73013E', 642: '#FF9079', 643: '#A79A74', 644: '#029BDB', 645: '#FF0169', 646: '#C7D2E7', 647: '#CA8869', 648: '#80FFCD', 649: '#BB1F69', 650: '#90B0AB', 651: '#7D74A9', 652: '#FCC7DB', 653: '#99375B', 654: '#00AB4D', 655: '#ABAED1', 656: '#BE9D91', 657: '#E6E5A7', 658: '#332C22', 659: '#DD587B', 660: '#F5FFF7', 661: '#5D3033', 662: '#6D3800', 663: '#FF0020', 664: '#B57BB3', 665: '#D7FFE6', 666: '#C535A9', 667: '#260009', 668: '#6A8781', 669: '#A8ABB4', 670: '#D45262', 671: '#794B61', 672: '#4621B2', 673: '#8DA4DB', 674: '#C7C890', 675: '#6FE9AD', 676: '#A243A7', 677: '#B2B081', 678: '#181B00', 679: '#286154', 680: '#4CA43B', 681: '#6A9573', 682: '#A8441D', 683: '#5C727B', 684: '#738671', 685: '#D0CFCB', 686: '#897B77', 687: '#1F3F22', 688: '#4145A7', 689: '#DA9894', 690: '#A1757A', 691: '#63243C', 692: '#ADAAFF', 693: '#00CDE2', 694: '#DDBC62', 695: '#698EB1', 696: '#208462', 697: '#00B7E0', 698: '#614A44', 699: '#9BBB57', 700: '#7A5C54', 701: '#857A50', 702: '#766B7E', 703: '#014833', 704: '#FF8347', 705: '#7A8EBA', 706: '#274740', 707: '#946444', 708: '#EBD8E6', 709: '#646241', 710: '#373917', 711: '#6AD450', 712: '#81817B', 713: '#D499E3', 714: '#979440', 715: '#011A12', 716: '#526554', 717: '#B5885C', 718: '#A499A5', 719: '#03AD89', 720: '#B3008B', 721: '#E3C4B5', 722: '#96531F', 723: '#867175', 724: '#74569E', 725: '#617D9F', 726: '#E70452', 727: '#067EAF', 728: '#A697B6', 729: '#B787A8', 730: '#9CFF93', 731: '#311D19', 732: '#3A9459', 733: '#6E746E', 734: '#B0C5AE', 735: '#84EDF7', 736: '#ED3488', 737: '#754C78', 738: '#384644', 739: '#C7847B', 740: '#00B6C5', 741: '#7FA670', 742: '#C1AF9E', 743: '#2A7FFF', 744: '#72A58C', 745: '#FFC07F', 746: '#9DEBDD', 747: '#D97C8E', 748: '#7E7C93', 749: '#62E674', 750: '#B5639E', 751: '#FFA861', 752: '#C2A580', 753: '#8D9C83', 754: '#B70546', 755: '#372B2E', 756: '#0098FF', 757: '#985975', 758: '#20204C', 759: '#FF6C60', 760: '#445083', 761: '#8502AA', 762: '#72361F', 763: '#9676A3', 764: '#484449', 765: '#CED6C2', 766: '#3B164A', 767: '#CCA763', 768: '#2C7F77', 769: '#02227B', 770: '#A37E6F', 771: '#CDE6DC', 772: '#CDFFFB', 773: '#BE811A', 774: '#F77183', 775: '#EDE6E2', 776: '#CDC6B4', 777: '#FFE09E', 778: '#3A7271', 779: '#FF7B59', 780: '#4E4E01', 781: '#4AC684', 782: '#8BC891', 783: '#BC8A96', 784: '#CF6353', 785: '#DCDE5C', 786: '#5EAADD', 787: '#F6A0AD', 788: '#E269AA', 789: '#A3DAE4', 790: '#436E83', 791: '#002E17', 792: '#ECFBFF', 793: '#A1C2B6', 794: '#50003F', 795: '#71695B', 796: '#67C4BB', 797: '#536EFF', 798: '#5D5A48', 799: '#890039', 800: '#969381', 801: '#371521', 802: '#5E4665', 803: '#AA62C3', 804: '#8D6F81', 805: '#2C6135', 806: '#410601', 807: '#564620', 808: '#E69034', 809: '#6DA6BD', 810: '#E58E56', 811: '#E3A68B', 812: '#48B176', 813: '#D27D67', 814: '#B5B268', 815: '#7F8427', 816: '#FF84E6', 817: '#435740', 818: '#EAE408', 819: '#F4F5FF', 820: '#325800', 821: '#4B6BA5', 822: '#ADCEFF', 823: '#9B8ACC', 824: '#885138', 825: '#5875C1', 826: '#7E7311', 827: '#FEA5CA', 828: '#9F8B5B', 829: '#A55B54', 830: '#89006A', 831: '#AF756F', 832: '#2A2000', 833: '#7499A1', 834: '#FFB550', 835: '#00011E', 836: '#D1511C', 837: '#688151', 838: '#BC908A', 839: '#78C8EB', 840: '#8502FF', 841: '#483D30', 842: '#C42221', 843: '#5EA7FF', 844: '#785715', 845: '#0CEA91', 846: '#FFFAED', 847: '#B3AF9D', 848: '#3E3D52', 849: '#5A9BC2', 850: '#9C2F90', 851: '#8D5700', 852: '#ADD79C', 853: '#00768B', 854: '#337D00', 855: '#C59700', 856: '#3156DC', 857: '#944575', 858: '#ECFFDC', 859: '#D24CB2', 860: '#97703C', 861: '#4C257F', 862: '#9E0366', 863: '#88FFEC', 864: '#B56481', 865: '#396D2B', 866: '#56735F', 867: '#988376', 868: '#9BB195', 869: '#A9795C', 870: '#E4C5D3', 871: '#9F4F67', 872: '#1E2B39', 873: '#664327', 874: '#AFCE78', 875: '#322EDF', 876: '#86B487', 877: '#C23000', 878: '#ABE86B', 879: '#96656D', 880: '#250E35', 881: '#A60019', 882: '#0080CF', 883: '#CAEFFF', 884: '#323F61', 885: '#A449DC', 886: '#6A9D3B', 887: '#FF5AE4', 888: '#636A01', 889: '#D16CDA', 890: '#736060', 891: '#FFBAAD', 892: '#D369B4', 893: '#FFDED6', 894: '#6C6D74', 895: '#927D5E', 896: '#845D70', 897: '#5B62C1', 898: '#2F4A36', 899: '#E45F35', 900: '#FF3B53', 901: '#AC84DD', 902: '#762988', 903: '#70EC98', 904: '#408543', 905: '#2C3533', 906: '#2E182D', 907: '#323925', 908: '#19181B', 909: '#2F2E2C', 910: '#023C32', 911: '#9B9EE2', 912: '#58AFAD', 913: '#5C424D', 914: '#7AC5A6', 915: '#685D75', 916: '#B9BCBD', 917: '#834357', 918: '#1A7B42', 919: '#2E57AA', 920: '#E55199', 921: '#316E47', 922: '#CD00C5', 923: '#6A004D', 924: '#7FBBEC', 925: '#F35691', 926: '#D7C54A', 927: '#62ACB7', 928: '#CBA1BC', 929: '#A28A9A', 930: '#6C3F3B', 931: '#FFE47D', 932: '#DCBAE3', 933: '#5F816D', 934: '#3A404A', 935: '#7DBF32', 936: '#E6ECDC', 937: '#852C19', 938: '#285366', 939: '#B8CB9C', 940: '#0E0D00', 941: '#4B5D56', 942: '#6B543F', 943: '#E27172', 944: '#0568EC', 945: '#2EB500', 946: '#D21656', 947: '#EFAFFF', 948: '#682021', 949: '#2D2011', 950: '#DA4CFF', 951: '#70968E', 952: '#FF7B7D', 953: '#4A1930', 954: '#E8C282', 955: '#E7DBBC', 956: '#A68486', 957: '#1F263C', 958: '#36574E', 959: '#52CE79', 960: '#ADAAA9', 961: '#8A9F45', 962: '#6542D2', 963: '#00FB8C', 964: '#5D697B', 965: '#CCD27F', 966: '#94A5A1', 967: '#790229', 968: '#E383E6', 969: '#7EA4C1', 970: '#4E4452', 971: '#4B2C00', 972: '#620B70', 973: '#314C1E', 974: '#874AA6', 975: '#E30091', 976: '#66460A', 977: '#EB9A8B', 978: '#EAC3A3', 979: '#98EAB3', 980: '#AB9180', 981: '#B8552F', 982: '#1A2B2F', 983: '#94DDC5', 984: '#9D8C76', 985: '#9C8333', 986: '#94A9C9', 987: '#392935', 988: '#8C675E', 989: '#CCE93A', 990: '#917100', 991: '#01400B', 992: '#449896', 993: '#1CA370', 994: '#E08DA7', 995: '#8B4A4E', 996: '#667776', 997: '#4692AD', 998: '#67BDA8', 999: '#69255C', 1000: '#D3BFFF', 1001: '#4A5132', 1002: '#7E9285', 1003: '#77733C', 1004: '#E7A0CC', 1005: '#51A288', 1006: '#2C656A', 1007: '#4D5C5E', 1008: '#C9403A', 1009: '#DDD7F3', 1010: '#005844', 1011: '#B4A200', 1012: '#488F69', 1013: '#858182', 1014: '#D4E9B9', 1015: '#3D7397', 1016: '#CAE8CE', 1017: '#D60034', 1018: '#AA6746', 1019: '#9E5585', 1020: '#BA6200'}
	#perform pca
	pca = PCA(n_components=components)
	pca.fit(array)
	explained = pca.explained_variance_ratio_
	#reduce dimensions
	X = pca.transform(array)
	#add colors
	label_color = [LABEL_COLOR_MAP[l] for l in labels]
	#probabilities sizes
	if type(probabilities) != int: size = 50 * probabilities.max(1)
	#plot
	plt.clf()
	plt.cla()
	plt.close()
	fig = plt.figure(1, figsize=(25, 22))
	if components == 2:
		#scatter plot
		if type(probabilities) != int:
			plt.scatter(X[:,0],X[:,1], c=label_color, s=size)
		else:
			plt.scatter(X[:,0],X[:,1], c=label_color)
		#label axes
		plt.xlabel('First principal component: %s' % explained[0], fontsize=20)
		plt.ylabel('Second principal component: %s' % explained[1], fontsize=20)
		#add title
		plt.title('2D Principal Component Analysis ' + title, fontsize=20)
		#add text
		for n, name in enumerate(names):
			plt.text(X[n,0], X[n,1], name, fontsize=10)
		#save or show
		plt.savefig(outfile)
	#plot 3d plot
	elif components == 3:
		#prepare figure
		#fig = plt.figure(1, figsize=(4, 3))
		plt.clf()
		ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
		plt.cla()
		#scatterplot
		if type(probabilities) != int:
			ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=label_color, cmap=plt.cm.spectral, edgecolor='k', s=size)
		else:
			ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=label_color, cmap=plt.cm.spectral, edgecolor='k')
		#label axes
		ax.set_xlabel('First principal component: %s' % explained[0])
		ax.set_ylabel('Second principal component: %s' % explained[1])
		ax.set_zlabel('Third principal component: %s' % explained[2])
		#add title
		ax.set_title('3D Principal Component Analysis ' + title)
		#add text
		for n, name in enumerate(names):
			ax.text(X[n,0], X[n,1], X[n,2], name, size=10)
		#save or show
		plt.show()
		#plt.savefig(outfile)
	else:
		raise Exception("Weird number of components")

def plot_heatmap(array, keys, patient_list, outfile, method, label, metric='euclidean', maximum=1.0):
	#clean up
	plt.clf()
	plt.cla()
	plt.close()
	sns.set(font_scale=0.25)
	#df = DataFrame(clean_array, columns=clean_keys, index=patient_list)
	df = DataFrame(array.T, index=keys, columns=patient_list)
	#legend_TN = [mpatches.Patch(color=c, label=l) for c,l in df[['tissue type','label']].drop_duplicates().values]
	if array.shape[0] != 1:
		g = sns.clustermap(df, col_cluster=True, row_cluster=False, metric=metric, cmap="Blues", vmax=maximum)
		g.fig.suptitle("Clustermap for %s, cluster %s" % (method, str(label+1)))
		g.ax_row_dendrogram.set_visible(False)
		g.ax_col_dendrogram.set_visible(False)
		#get row order and col order
		column_order = [df.columns[i] for i in g.dendrogram_col.reordered_ind]
		#g.cax.set_visible(False)
		plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
	else:
		column_order = df.columns
		if "binary" not in method:
			#g = sns.heatmap(df, cbar=False)
			g = sns.heatmap(df, cmap="Blues", vmax=maximum)
		else:
			#g = sns.heatmap(df, cbar=False, vmin=0, vmax=1)
			g = sns.heatmap(df, vmin=0, vmax=1, cmap="Blues")
		g.set_title("Clustermap for %s, cluster %s" % (method, str(label+1)))
		plt.yticks(rotation=0)
	#g.set(xlabel='my x label', ylabel='my y label')
	heatmap_outfile = outfile.rsplit(".",1)[0] + "_" + str(label+1) + "." + outfile.rsplit(".",1)[1]
	#plt.savefig(heatmap_outfile)
	#return df with new column order
	return df.reindex(columns=column_order)

def plot_clusters(labels, array, keys, patient_list, outfile, method, metric='euclidean'):
	"""
	Plot a heatmap for each cluster
	"""
	clean_vector_array = np.array(array)
	clean_vector_keys = list(keys)
	#define maximum
	if "binary" in method: maximum = 1.0
	else: maximum = max(list(itertools.chain.from_iterable(array.tolist())))
	flag = 1
	#remove 0s
	while flag == 1:
		flag = 0
		for n, column in enumerate(clean_vector_array.T):
			if not np.any(column):
				clean_vector_array = np.delete(clean_vector_array,n,1)
				del clean_vector_keys[n]
				flag = 1
				break
	#plot a heatmap/clustergram for each cluster
	label_patients_dict = {}
	dfs = {}
	#get proper labels
	the_labels = list(set(labels))
	#iterate through labels
	for label in the_labels:
		label_patients_dict[label] = []
		#get relevant indices
		indices = []
		label_patients = []
		for n,i in enumerate(labels):
			if i == label:
				indices.append(n)
				label_patients_dict[label].append(patient_list[n])
		#relevant patients
		label_array = np.array(clean_vector_array[indices])
		label_keys = list(clean_vector_keys)
		#plot
		dfs[label] = plot_heatmap(label_array, label_keys, label_patients_dict[label], outfile, method, label, metric, maximum)
	#colors
	colors_limited = ["Blues", "Oranges", "Greens", "Purples", "Reds", "Greys"]
	colors = ["Blues", "Oranges", "Greens", "Purples", "Reds", "Greys", 'YlGn', 'BrBG', 'OrRd','BuGn', 'BuPu', 'CMRmap', 'GnBu', 'PRGn', 'PiYG', 'PuBu',
			  'PuBuGn', 'PuOr', 'PuRd', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Spectral', 'Wistia', 'YlGnBu', 'YlOrBr',
			  'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cool', 'coolwarm', 'copper', 'cubehelix', 'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow',
			  'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'icefire', 'inferno', 'jet', 'magma', 'mako', 'nipy_spectral', 'ocean', 'pink', 'plasma', 'rainbow',
			  'rocket', 'seismic', 'spectral', 'spring', 'summer', 'terrain', 'viridis', 'vlag', 'winter']
	#plot all together in one heatmap
	#clear
	plt.clf()
	plt.cla()
	plt.close()
	#initialize
	#limit labels
	limited_labels = [label for label in the_labels if len(label_patients_dict[label]) > 1]
	#sort labels
	limited_labels = sorted(limited_labels, key=lambda x:len(label_patients_dict[x]),reverse=True)
	#limit only to top clusters
	MAX_CLUSTERS = 10
	if len(limited_labels) > MAX_CLUSTERS: limited_labels = limited_labels[:MAX_CLUSTERS]
	#create figure
	fig = plt.figure(1, figsize=(6.3, 6.3), dpi=300)
	#create subplots with proportion according to number of patients in each label
	gs = gridspec.GridSpec(2, len(limited_labels), height_ratios=[1,40], width_ratios=[len(label_patients_dict[label]) for label in limited_labels], wspace = 0.01, hspace=0.03)
	gs_cbar = gridspec.GridSpec(2, len(limited_labels), height_ratios=[1,40], wspace = 0.03, hspace=0.03)
	#initialize axes
	axs = {}
	cbar_axs = {}
	for n in range(2*len(limited_labels)):
		if n < len(limited_labels):
			ax = plt.subplot(gs_cbar[n])
			cbar_axs[limited_labels[n]] = ax
			ax.get_xaxis().set_visible(False)
			ax.get_yaxis().set_visible(False)
		else:
			ax = plt.subplot(gs[n])
			axs[limited_labels[n - len(limited_labels)]] = ax
	#iterate through labels and plot
	for n, label in enumerate(limited_labels):
		#get color
		color = n%len(colors)
		#get yticklabels
		if n == 0: yticklabels = True #only first label gest yticks
		else: yticklabels = False
		#get ratio for color bar
		#ratio = len(the_labels)*float(2)/ len(label_patients_dict[label])
		ratio = float(2)/ len(label_patients_dict[label])
		#width = str(ratio*100) + "%"
		width = "90%"
		#plot heatmap!
		cax = inset_axes(cbar_axs[label], width=width, height="50%", loc=10)
		g = sns.heatmap(dfs[label], cmap=colors[color], cbar=True, cbar_ax=cax, cbar_kws = {"orientation": "horizontal", "ticks":[0.0,maximum]}, vmax=maximum, ax=axs[label], yticklabels=yticklabels, xticklabels=True)
		#fix color bar ticks
		cax.tick_params(width=0.5, length=1, labelsize=4, direction="out", labeltop='on', labelbottom='off', top='on', bottom='off')
		cax.set_xticklabels([str(round(float(i), 2)) for i in [0.0, maximum]])
		#change tick sizes x axis
		xticksize = min(float(200) / len(patient_list), 4)
		xticksize = max(xticksize, 0.5)
		g.set_xticklabels(g.get_xmajorticklabels(), fontsize = xticksize)
		#y axis
		#g.set_yticklabels(g.get_ymajorticklabels(), fontsize = 1.25)
		yticksize = min(float(300) / len(clean_vector_keys), 5)
		g.set_yticklabels(g.get_ymajorticklabels(), fontsize = yticksize)
		#annotate
		plt.setp(g.get_xticklabels(), rotation=90)
		plt.setp(g.get_yticklabels(), rotation=0)
		
	plt.subplots_adjust(left=0.25, right=0.95, bottom=0.12, top=0.98)
	
	#plotfile
	plotfile = outfile.rsplit(".",1)[0] + "_" + "ALL" + "." + outfile.rsplit(".",1)[1]
	#save
	plt.savefig(plotfile, transparent=True)

if __name__ == "__main__":
		exit()
