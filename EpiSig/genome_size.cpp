#include<cstdio>
#include<cstdlib>
#include<cstring>
#include "genome_size.h"
#include "utils.h"

void GS_init( GS& gs, string& name )
{
	GS_init( gs, name.c_str() );
}

void GS_init( GS& gs, const char* name )
{
	if( strcmp(name,"hg18") == 0 )
	{
		printf( "Using pre-defined hg18 genome assembly information\n" );
		gs.g_cnt = 24;

		gs.id_chr_map.push_back("chr1");
		gs.id_chr_map.push_back("chr2");
		gs.id_chr_map.push_back("chr3");
		gs.id_chr_map.push_back("chr4");
		gs.id_chr_map.push_back("chr5");
		gs.id_chr_map.push_back("chr6");
		gs.id_chr_map.push_back("chr7");
		gs.id_chr_map.push_back("chr8");
		gs.id_chr_map.push_back("chr9");
		gs.id_chr_map.push_back("chr10");
		gs.id_chr_map.push_back("chr11");
		gs.id_chr_map.push_back("chr12");
		gs.id_chr_map.push_back("chr13");
		gs.id_chr_map.push_back("chr14");
		gs.id_chr_map.push_back("chr15");
		gs.id_chr_map.push_back("chr16");
		gs.id_chr_map.push_back("chr17");
		gs.id_chr_map.push_back("chr18");
		gs.id_chr_map.push_back("chr19");
		gs.id_chr_map.push_back("chr20");
		gs.id_chr_map.push_back("chr21");
		gs.id_chr_map.push_back("chr22");
		gs.id_chr_map.push_back("chrX");
		gs.id_chr_map.push_back("chrY");

		gs.g_size.push_back(247249719); // chr1
		gs.g_size.push_back(242951149); // chr2
		gs.g_size.push_back(199501827); // chr3
		gs.g_size.push_back(191273063); // chr4
		gs.g_size.push_back(180857866); // chr5
		gs.g_size.push_back(170899992); // chr6
		gs.g_size.push_back(158821424); // chr7
		gs.g_size.push_back(146274826); // chr8
		gs.g_size.push_back(140273252); // chr9
		gs.g_size.push_back(135374737); // chr10
		gs.g_size.push_back(134452384); // chr11
		gs.g_size.push_back(132349534); // chr12
		gs.g_size.push_back(114142980); // chr13
		gs.g_size.push_back(106368585); // chr14
		gs.g_size.push_back(100338915); // chr15
		gs.g_size.push_back(88827254); // chr16
		gs.g_size.push_back(78774742); // chr17
		gs.g_size.push_back(76117153); // chr18
		gs.g_size.push_back(63811651); // chr19
		gs.g_size.push_back(62435964); // chr20
		gs.g_size.push_back(46944323); // chr21
		gs.g_size.push_back(49691432); // chr22
		gs.g_size.push_back(154913754); // chrX
		gs.g_size.push_back(57772954); // chrY

		gs.chr_id_map["chr1"] = 0;
		gs.chr_id_map["chr2"] = 1;
		gs.chr_id_map["chr3"] = 2;
		gs.chr_id_map["chr4"] = 3;
		gs.chr_id_map["chr5"] = 4;
		gs.chr_id_map["chr6"] = 5;
		gs.chr_id_map["chr7"] = 6;
		gs.chr_id_map["chr8"] = 7;
		gs.chr_id_map["chr9"] = 8;
		gs.chr_id_map["chr10"] = 9;
		gs.chr_id_map["chr11"] = 10;
		gs.chr_id_map["chr12"] = 11;
		gs.chr_id_map["chr13"] = 12;
		gs.chr_id_map["chr14"] = 13;
		gs.chr_id_map["chr15"] = 14;
		gs.chr_id_map["chr16"] = 15;
		gs.chr_id_map["chr17"] = 16;
		gs.chr_id_map["chr18"] = 17;
		gs.chr_id_map["chr19"] = 18;
		gs.chr_id_map["chr20"] = 19;
		gs.chr_id_map["chr21"] = 20;
		gs.chr_id_map["chr22"] = 21;
		gs.chr_id_map["chrX"] = 22;
		gs.chr_id_map["chrY"] = 23;
	}
	else if( strcmp(name,"hg19") == 0 )
	{
		printf( "Using pre-defined hg19 genome assembly information\n" );
		gs.g_cnt = 24;

		gs.id_chr_map.push_back("chr1");
		gs.id_chr_map.push_back("chr2");
		gs.id_chr_map.push_back("chr3");
		gs.id_chr_map.push_back("chr4");
		gs.id_chr_map.push_back("chr5");
		gs.id_chr_map.push_back("chr6");
		gs.id_chr_map.push_back("chr7");
		gs.id_chr_map.push_back("chr8");
		gs.id_chr_map.push_back("chr9");
		gs.id_chr_map.push_back("chr10");
		gs.id_chr_map.push_back("chr11");
		gs.id_chr_map.push_back("chr12");
		gs.id_chr_map.push_back("chr13");
		gs.id_chr_map.push_back("chr14");
		gs.id_chr_map.push_back("chr15");
		gs.id_chr_map.push_back("chr16");
		gs.id_chr_map.push_back("chr17");
		gs.id_chr_map.push_back("chr18");
		gs.id_chr_map.push_back("chr19");
		gs.id_chr_map.push_back("chr20");
		gs.id_chr_map.push_back("chr21");
		gs.id_chr_map.push_back("chr22");
		gs.id_chr_map.push_back("chrX");
		gs.id_chr_map.push_back("chrY");

		gs.g_size.push_back(249250621); // chr1
		gs.g_size.push_back(243199373); // chr2
		gs.g_size.push_back(198022430); // chr3
		gs.g_size.push_back(191154276); // chr4
		gs.g_size.push_back(180915260); // chr5
		gs.g_size.push_back(171115067); // chr6
		gs.g_size.push_back(159138663); // chr7
		gs.g_size.push_back(146364022); // chr8
		gs.g_size.push_back(141213431); // chr9
		gs.g_size.push_back(135534747); // chr10
		gs.g_size.push_back(135006516); // chr11
		gs.g_size.push_back(133851895); // chr12
		gs.g_size.push_back(115169878); // chr13
		gs.g_size.push_back(107349540); // chr14
		gs.g_size.push_back(102531392); // chr15
		gs.g_size.push_back(90354753); // chr16
		gs.g_size.push_back(81195210); // chr17
		gs.g_size.push_back(78077248); // chr18
		gs.g_size.push_back(59128983); // chr19
		gs.g_size.push_back(63025520); // chr20
		gs.g_size.push_back(48129895); // chr21
		gs.g_size.push_back(51304566); // chr22
		gs.g_size.push_back(155270560); // chrX
		gs.g_size.push_back(59373566); // chrY

		gs.chr_id_map["chr1"] = 0;
		gs.chr_id_map["chr2"] = 1;
		gs.chr_id_map["chr3"] = 2;
		gs.chr_id_map["chr4"] = 3;
		gs.chr_id_map["chr5"] = 4;
		gs.chr_id_map["chr6"] = 5;
		gs.chr_id_map["chr7"] = 6;
		gs.chr_id_map["chr8"] = 7;
		gs.chr_id_map["chr9"] = 8;
		gs.chr_id_map["chr10"] = 9;
		gs.chr_id_map["chr11"] = 10;
		gs.chr_id_map["chr12"] = 11;
		gs.chr_id_map["chr13"] = 12;
		gs.chr_id_map["chr14"] = 13;
		gs.chr_id_map["chr15"] = 14;
		gs.chr_id_map["chr16"] = 15;
		gs.chr_id_map["chr17"] = 16;
		gs.chr_id_map["chr18"] = 17;
		gs.chr_id_map["chr19"] = 18;
		gs.chr_id_map["chr20"] = 19;
		gs.chr_id_map["chr21"] = 20;
		gs.chr_id_map["chr22"] = 21;
		gs.chr_id_map["chrX"] = 22;
		gs.chr_id_map["chrY"] = 23;
	}
	else if( strcmp(name,"hg38") == 0 )
	{
		printf( "Using pre-defined hg38 genome assembly information\n" );
		gs.g_cnt = 24;

		gs.id_chr_map.push_back("chr1");
		gs.id_chr_map.push_back("chr2");
		gs.id_chr_map.push_back("chr3");
		gs.id_chr_map.push_back("chr4");
		gs.id_chr_map.push_back("chr5");
		gs.id_chr_map.push_back("chr6");
		gs.id_chr_map.push_back("chr7");
		gs.id_chr_map.push_back("chr8");
		gs.id_chr_map.push_back("chr9");
		gs.id_chr_map.push_back("chr10");
		gs.id_chr_map.push_back("chr11");
		gs.id_chr_map.push_back("chr12");
		gs.id_chr_map.push_back("chr13");
		gs.id_chr_map.push_back("chr14");
		gs.id_chr_map.push_back("chr15");
		gs.id_chr_map.push_back("chr16");
		gs.id_chr_map.push_back("chr17");
		gs.id_chr_map.push_back("chr18");
		gs.id_chr_map.push_back("chr19");
		gs.id_chr_map.push_back("chr20");
		gs.id_chr_map.push_back("chr21");
		gs.id_chr_map.push_back("chr22");
		gs.id_chr_map.push_back("chrX");
		gs.id_chr_map.push_back("chrY");

		gs.g_size.push_back(248956422); // chr1
		gs.g_size.push_back(242193529); // chr2
		gs.g_size.push_back(198295559); // chr3
		gs.g_size.push_back(190214555); // chr4
		gs.g_size.push_back(181538259); // chr5
		gs.g_size.push_back(170805979); // chr6
		gs.g_size.push_back(159345973); // chr7
		gs.g_size.push_back(145138636); // chr8
		gs.g_size.push_back(138394717); // chr9
		gs.g_size.push_back(133797422); // chr10
		gs.g_size.push_back(135086622); // chr11
		gs.g_size.push_back(133275309); // chr12
		gs.g_size.push_back(114364328); // chr13
		gs.g_size.push_back(107043718); // chr14
		gs.g_size.push_back(101991189); // chr15
		gs.g_size.push_back(90338345); // chr16
		gs.g_size.push_back(83257441); // chr17
		gs.g_size.push_back(80373285); // chr18
		gs.g_size.push_back(58617616); // chr19
		gs.g_size.push_back(64444167); // chr20
		gs.g_size.push_back(46709983); // chr21
		gs.g_size.push_back(50818468); // chr22
		gs.g_size.push_back(156040895); // chrX
		gs.g_size.push_back(57227415); // chrY

		gs.chr_id_map["chr1"] = 0;
		gs.chr_id_map["chr2"] = 1;
		gs.chr_id_map["chr3"] = 2;
		gs.chr_id_map["chr4"] = 3;
		gs.chr_id_map["chr5"] = 4;
		gs.chr_id_map["chr6"] = 5;
		gs.chr_id_map["chr7"] = 6;
		gs.chr_id_map["chr8"] = 7;
		gs.chr_id_map["chr9"] = 8;
		gs.chr_id_map["chr10"] = 9;
		gs.chr_id_map["chr11"] = 10;
		gs.chr_id_map["chr12"] = 11;
		gs.chr_id_map["chr13"] = 12;
		gs.chr_id_map["chr14"] = 13;
		gs.chr_id_map["chr15"] = 14;
		gs.chr_id_map["chr16"] = 15;
		gs.chr_id_map["chr17"] = 16;
		gs.chr_id_map["chr18"] = 17;
		gs.chr_id_map["chr19"] = 18;
		gs.chr_id_map["chr20"] = 19;
		gs.chr_id_map["chr21"] = 20;
		gs.chr_id_map["chr22"] = 21;
		gs.chr_id_map["chrX"] = 22;
		gs.chr_id_map["chrY"] = 23;
	}
	else if( strcmp(name,"mm9") == 0 )
	{
		printf( "Using pre-defined mm9 genome assembly information\n" );
		gs.g_cnt = 21;

		gs.id_chr_map.push_back("chr1");
		gs.id_chr_map.push_back("chr2");
		gs.id_chr_map.push_back("chr3");
		gs.id_chr_map.push_back("chr4");
		gs.id_chr_map.push_back("chr5");
		gs.id_chr_map.push_back("chr6");
		gs.id_chr_map.push_back("chr7");
		gs.id_chr_map.push_back("chr8");
		gs.id_chr_map.push_back("chr9");
		gs.id_chr_map.push_back("chr10");
		gs.id_chr_map.push_back("chr11");
		gs.id_chr_map.push_back("chr12");
		gs.id_chr_map.push_back("chr13");
		gs.id_chr_map.push_back("chr14");
		gs.id_chr_map.push_back("chr15");
		gs.id_chr_map.push_back("chr16");
		gs.id_chr_map.push_back("chr17");
		gs.id_chr_map.push_back("chr18");
		gs.id_chr_map.push_back("chr19");
		gs.id_chr_map.push_back("chrX");
		gs.id_chr_map.push_back("chrY");

		gs.g_size.push_back(197195432); // chr1
		gs.g_size.push_back(181748087); // chr2
		gs.g_size.push_back(159599783); // chr3
		gs.g_size.push_back(155630120); // chr4
		gs.g_size.push_back(152537259); // chr5
		gs.g_size.push_back(149517037); // chr6
		gs.g_size.push_back(152524553); // chr7
		gs.g_size.push_back(131738871); // chr8
		gs.g_size.push_back(124076172); // chr9
		gs.g_size.push_back(129993255); // chr10
		gs.g_size.push_back(121843856); // chr11
		gs.g_size.push_back(121257530); // chr12
		gs.g_size.push_back(120284312); // chr13
		gs.g_size.push_back(125194864); // chr14
		gs.g_size.push_back(103494974); // chr15
		gs.g_size.push_back(98319150); // chr16
		gs.g_size.push_back(95272651); // chr17
		gs.g_size.push_back(90772031); // chr18
		gs.g_size.push_back(61342430); // chr19
		gs.g_size.push_back(166650296); // chrX
		gs.g_size.push_back(15902555); // chrY

		gs.chr_id_map["chr1"] = 0;
		gs.chr_id_map["chr2"] = 1;
		gs.chr_id_map["chr3"] = 2;
		gs.chr_id_map["chr4"] = 3;
		gs.chr_id_map["chr5"] = 4;
		gs.chr_id_map["chr6"] = 5;
		gs.chr_id_map["chr7"] = 6;
		gs.chr_id_map["chr8"] = 7;
		gs.chr_id_map["chr9"] = 8;
		gs.chr_id_map["chr10"] = 9;
		gs.chr_id_map["chr11"] = 10;
		gs.chr_id_map["chr12"] = 11;
		gs.chr_id_map["chr13"] = 12;
		gs.chr_id_map["chr14"] = 13;
		gs.chr_id_map["chr15"] = 14;
		gs.chr_id_map["chr16"] = 15;
		gs.chr_id_map["chr17"] = 16;
		gs.chr_id_map["chr18"] = 17;
		gs.chr_id_map["chr19"] = 18;
		gs.chr_id_map["chrX"] = 19;
		gs.chr_id_map["chrY"] = 20;
	}
	else if( strcmp(name,"mm10") == 0 )
	{
		printf( "Using pre-defined mm10 genome assembly information\n" );
		gs.g_cnt = 21;

		gs.id_chr_map.push_back("chr1");
		gs.id_chr_map.push_back("chr2");
		gs.id_chr_map.push_back("chr3");
		gs.id_chr_map.push_back("chr4");
		gs.id_chr_map.push_back("chr5");
		gs.id_chr_map.push_back("chr6");
		gs.id_chr_map.push_back("chr7");
		gs.id_chr_map.push_back("chr8");
		gs.id_chr_map.push_back("chr9");
		gs.id_chr_map.push_back("chr10");
		gs.id_chr_map.push_back("chr11");
		gs.id_chr_map.push_back("chr12");
		gs.id_chr_map.push_back("chr13");
		gs.id_chr_map.push_back("chr14");
		gs.id_chr_map.push_back("chr15");
		gs.id_chr_map.push_back("chr16");
		gs.id_chr_map.push_back("chr17");
		gs.id_chr_map.push_back("chr18");
		gs.id_chr_map.push_back("chr19");
		gs.id_chr_map.push_back("chrX");
		gs.id_chr_map.push_back("chrY");

		gs.g_size.push_back(195471971); // chr1
		gs.g_size.push_back(182113224); // chr2
		gs.g_size.push_back(160039680); // chr3
		gs.g_size.push_back(156508116); // chr4
		gs.g_size.push_back(151834684); // chr5
		gs.g_size.push_back(149736546); // chr6
		gs.g_size.push_back(145441459); // chr7
		gs.g_size.push_back(129401213); // chr8
		gs.g_size.push_back(124595110); // chr9
		gs.g_size.push_back(130694993); // chr10
		gs.g_size.push_back(122082543); // chr11
		gs.g_size.push_back(120129022); // chr12
		gs.g_size.push_back(120421639); // chr13
		gs.g_size.push_back(124902244); // chr14
		gs.g_size.push_back(104043685); // chr15
		gs.g_size.push_back(98207768); // chr16
		gs.g_size.push_back(94987271); // chr17
		gs.g_size.push_back(90702639); // chr18
		gs.g_size.push_back(61431566); // chr19
		gs.g_size.push_back(171031299); // chrX
		gs.g_size.push_back(91744698); // chrY

		gs.chr_id_map["chr1"] = 0;
		gs.chr_id_map["chr2"] = 1;
		gs.chr_id_map["chr3"] = 2;
		gs.chr_id_map["chr4"] = 3;
		gs.chr_id_map["chr5"] = 4;
		gs.chr_id_map["chr6"] = 5;
		gs.chr_id_map["chr7"] = 6;
		gs.chr_id_map["chr8"] = 7;
		gs.chr_id_map["chr9"] = 8;
		gs.chr_id_map["chr10"] = 9;
		gs.chr_id_map["chr11"] = 10;
		gs.chr_id_map["chr12"] = 11;
		gs.chr_id_map["chr13"] = 12;
		gs.chr_id_map["chr14"] = 13;
		gs.chr_id_map["chr15"] = 14;
		gs.chr_id_map["chr16"] = 15;
		gs.chr_id_map["chr17"] = 16;
		gs.chr_id_map["chr18"] = 17;
		gs.chr_id_map["chr19"] = 18;
		gs.chr_id_map["chrX"] = 19;
		gs.chr_id_map["chrY"] = 20;
	}
	else
	{
		printf( "Using user defined genome assembly information from %s\n", name );
		GS_init_preset( gs, name );
	}
}

void GS_init_preset( GS& gs, const char* name )
{
	int cnt = 2;
	FILE* fp;
	char line[100];
	string fields[2];
	
	if( (fp=fopen(name,"r")) == NULL )
	{
		printf( "Failed to open file: %s\n", name );
		exit(1);
	}

	gs.g_cnt = 0;
	while( fgets(line,99,fp) != NULL )
	{
		spliter(line,'\t',fields,cnt);

		gs.g_size.push_back(atol(fields[1].c_str()));
		gs.id_chr_map.push_back(fields[0]);
		gs.chr_id_map[fields[0]] = gs.g_cnt;

		gs.g_cnt++;
	}

	fclose(fp);
}
