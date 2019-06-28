#ifndef GAURD_DECS_1
#define GAURD_DECS_1

// marked2

std::vector<double> linspace(const double start,const double stop,const size_t N){
	using namespace std;
	if(stop<start)
		throw(domain_error("ERROR in linspace"));
	
	const double step( (stop-start)/(N-1) );
	
	vector<double> val(N);
	val[0]=start;
	for(size_t i(1);i<N;++i){
		val[i]=val[i-1]+step;
	}
	return val;
}


std::vector<double> range(const size_t N){
	std::vector<double> A(N,0);
	for(size_t i(0);i<N;++i)
		A[i]=i;
	return A;
}


size_t ALT_count_to_combos(const size_t ALT_count){
	return (ALT_count+1)*(ALT_count+2)/2;
}


size_t combos_to_ALT_count(const size_t C){
	return (pow((C-1)*8+9,0.5)-3)/2;
}



std::vector<std::string> split(const std::string& line){
    std::vector<std::string> split_string;
    std::string::const_iterator i(line.begin()),j(line.begin());
    while(i!=line.end() && j!=line.end()){
        while(isspace(*i) && i!=line.end())
            ++i;
        j=i;
        while(!isspace(*j) && j!=line.end())
            ++j;
        std::string word(i,j);
        if(word.length()>0)
            split_string.push_back(word);
        i=j;
    }
    return split_string;
}


std::string shorten_word(const std::string& word){
	std::string::const_iterator i(word.begin()),j;
	while(isspace(*i) && i!=word.end())
		++i;
	j=i;
	while(!isspace(*j) && j!=word.end())
		++j;
	return std::string(i,j);
}


std::vector<std::string> split(const std::string& line,const char sep){
	std::vector<std::string> split_string;
	std::string::const_iterator i(line.begin()),j(line.begin());
	while(i!=line.end() && j!=line.end()){
		while(*j!=sep && j!=line.end())
			++j;
		std::string full_word(std::string(i,j));
		//std::cout<<full_word.size()<<'\n';
		std::string word(shorten_word(full_word));
		split_string.push_back(word);
		if(j!=line.end())
			i=j+1;
		else
			i=j;
		j=i;
	}
	if(line[line.size()-1]==sep)
		split_string.push_back("");
	return split_string;
}






void col_extractor(const std::string vcf_header,size_t& CHROM_col,size_t& POS_col,size_t& ID_col,size_t& REF_col,size_t& ALT_col,size_t& QUAL_col,size_t& FILTER_col,size_t& INFO_col,size_t& FORMAT_col){
	std::vector<std::string> split_line(split(vcf_header));
	CHROM_col=POS_col=ID_col=REF_col=ALT_col=QUAL_col=FILTER_col=INFO_col=FORMAT_col=1000;
	for(size_t i(0);i<split_line.size();++i){
		std::string F(split_line[i]);
		if(F=="#CHROM")
			CHROM_col=i;
		if(F=="POS")
			POS_col=i;
		if(F=="ID")
			ID_col=i;
		if(F=="REF")
			REF_col=i;
		if(F=="ALT")
			ALT_col=i;
		if(F=="QUAL")
			QUAL_col=i;
		if(F=="FILTER")
			FILTER_col=i;
		if(F=="INFO")
			INFO_col=i;
		if(F=="FORMAT")
			FORMAT_col=i;
	}
}





template<class T>
std::vector<T> LOG_1D(const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=std::log(v[i]);
	return ans;
}


template<class T>
std::vector<std::vector<T> > LOG_2D(const std::vector<std::vector<T> >& v){
	size_t d1(v.size()),d2(v[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=log(v[i][j]);
		}
	}
	return ans;
}


template<class T>
std::vector<T> EXP_1D(const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=std::exp(v[i]);
	return ans;
}

template<class T>
std::vector<std::vector<T> > EXP_2D(const std::vector<std::vector<T> >& v){
	size_t d1(v.size()),d2(v[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=std::exp(v[i][j]);
		}
	}
	return ans;
}

// no need for row extract since = operator copies the entire row vector
template<class T>
std::vector<T> column_extract(const std::vector<std::vector<T> >& P,const size_t C){
	if(C>P[0].size()-1)
		throw(std::domain_error("error in column_extract"));
	const size_t states(P.size());
	std::vector<T> col(states,P[0][0]);
	for(size_t i(0);i<states;++i)
		col[i]=P[i][C];
	return col;
}

template<class T>
void column_insert(std::vector<std::vector<T> >& P,const std::vector<T>& Col,const size_t C){
	if(P.size()!=Col.size())
		throw(std::domain_error("error in column_insert"));
	const size_t states(P.size());
	for(size_t i(0);i<states;++i)
		P[i][C]=Col[i];
}



template<class T>
T sum(const std::vector<T>& v){
	T ans(0);
	for(size_t i(0);i<v.size();++i)
		ans+=v[i];
	return ans;
}

template<class T>
T sum_2D(const std::vector<std::vector<T> >& v){
	T ans(0);
	for(size_t i(0);i<v.size();++i)
		ans+=sum(v[i]);
	return ans;
}


template<class T>
T max(const std::vector<T>& v){
	T ans(v[0]);
	for(size_t i(1);i<v.size();++i)
		if(v[i]>ans)
			ans=v[i];
	return ans;
}

template<class T>
T max_2D(const std::vector<std::vector<T > >& v){
	using namespace std;
	T ans(max(v[0]));
	for(size_t i(1);i<v.size();++i){
		const T M( max(v[i]) );
		if(M>ans)
			ans=M;
	}
	return ans;
}

template<class T>
T max_abs(const std::vector<T>& v){
	T ans(std::abs(v[0]));
	for(size_t i(1);i<v.size();++i)
		if(std::abs(v[i])>ans)
			ans=std::abs(v[i]);
	return ans;
}


template<class T>
size_t arg_max(const std::vector<T>& v){
	size_t ans(0);
	for(size_t i(1);i<v.size();++i)
		if(v[i]>v[ans])
			ans=i;
	return ans;
}


template<class T>
size_t arg_max_abs(const std::vector<T>& v){
	size_t ans(0);
	for(size_t i(1);i<v.size();++i)
		if(std::abs(v[i])>std::abs(v[ans]))
			ans=i;
	return ans;
}



template<class T>
T min(const std::vector<T>& v){
	T ans(v[0]);
	for(size_t i(1);i<v.size();++i)
		if(v[i]<ans)
			ans=v[i];
	return ans;
}

template<class T>
T min_2D(const std::vector<std::vector<T> >& v){
	using namespace std;
	T ans(min(v[0]));
	for(size_t i(1);i<v.size();++i)
		if( min(v[i])<ans )
			ans=min(v[i]);
	return ans;
}



template<class T>
std::vector<T> dot(const std::vector<std::vector<T> >& v1,const std::vector<T>& v2){
	if(v1[0].size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v1.size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=sum(v1[i]*v2);
	return ans;
}

template<class T>
std::vector<T> dot(const std::vector<T>& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v2[0].size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=sum(column_extract(v2,i)*v1);
	return ans;
}


template<class T>
std::vector<std::vector<T> > dot(const std::vector<std::vector<T> >& v1,const std::vector<std::vector<T> >& v2){
	if(v1[0].size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	std::vector<std::vector<T> > ans(v1.size(),std::vector<T>(v2[0].size()));
	for(size_t i(0);i<ans.size();++i){
		for(size_t j(0);j<ans[0].size();++j){
			T summ(0);
			for(size_t k(0);k<v1[0].size();++k){
				summ+=v1[i][k]*v2[k][j];
			}
			ans[i][j]=summ;
		}
	}
	return ans;
}



template<class T>
std::vector<T> max_sum_dot(const std::vector<std::vector<T> >& v1,const std::vector<T>& v2){
	if(v1[0].size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v1.size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=max(v1[i]+v2);
	return ans;
}

template<class T>
std::vector<T> max_sum_dot(const std::vector<T>& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v2[0].size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=max(column_extract(v2,i)+v1);
	return ans;
}




template<class T>
std::vector<T> operator*(const double S,const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<v.size();++i)
		ans[i]=v[i]*S;
	return ans;
}



template<class T>
std::vector<T> operator*(const std::vector<T>& v,const double S){
	return S*v;
}

template<class T>
std::vector<T> operator+(const double S,const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<v.size();++i)
		ans[i]=v[i]+S;
	return ans;
}

template<class T>
std::vector<T> operator+(const std::vector<T>& v,const double S){
	return S+v;
}


template<class T>
std::vector<T> operator*(const std::vector<T>& v1,const std::vector<T>& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("sizes do not match in multiplication"));
	
	size_t l(v1.size());
	std::vector<T> ans(l,v1[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=v1[i]*v2[i];
	return ans;
}


template<class T>
std::vector<std::vector<T> > operator*(const std::vector<std::vector<T> >& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size() || v1[0].size()!=v2[0].size())
		throw(std::domain_error("sizes do not match in multiplication"));
	
	size_t d1(v1.size()),d2(v1[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v1[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=v1[i][j]*v2[i][j];
		}
	}
	return ans;
}


template<class T>
std::vector<T> operator+(const std::vector<T>& v1,const std::vector<T>& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("sizes do not match in addition"));
	
	size_t l(v1.size());
	std::vector<T> ans(l,v1[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=v1[i]+v2[i];
	return ans;
}


template<class T>
std::vector<std::vector<T> > operator+(const std::vector<std::vector<T> >& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size() || v1[0].size()!=v2[0].size())
		throw(std::domain_error("sizes do not match in addition"));
	
	size_t d1(v1.size()),d2(v1[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v1[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=v1[i][j]+v2[i][j];
		}
	}
	return ans;
}

template<class T1,class T2>
std::ostream& operator<<(std::ostream& out,const std::map<T1,T2>& V){
	using namespace std;
	typename map<T1,T2>::const_iterator I( V.begin() );
	for(;I!=V.end();++I){
		out<< I->first <<"\n";
		out<<"-----\n";
		out<< I->second <<"\n";
	}
	return out;
}



template<class T>
std::ostream& operator<<(std::ostream& out,const std::vector<T>& v){
	if(v.size()>0){
		for(size_t i(0);i<v.size();++i)
			out<<v[i]<<':';
	}
	return out;
}

template<class T>
std::ostream& operator<<(std::ostream& out,const std::list<T>& v){
	typename std::list<T>::const_iterator I(v.begin());
	while(I!=v.end()){
		std::cout<<(*I)<<":";
		++I;
	}
	return out;
}


template<class T>
std::ostream& operator<<(std::ostream& out,const std::vector<std::vector<T> >& v){
	if(v.size()>0){
		for(size_t i(0);i<v.size();++i)
			out<<v[i]<<'\n';
	}
	return out;
}





void col_extractor_for_FORMAT(const std::string& FORMAT,size_t& GT_col,size_t& AD_col,size_t& DP_col,size_t& GQ_col,size_t& PGT_col,size_t& PID_col,size_t& PL_col){
	GT_col=AD_col=DP_col=GQ_col=PGT_col=PID_col=PL_col=1000;
	std::vector<std::string> split_FORMAT(split(FORMAT,':'));
	for(size_t i(0);i<split_FORMAT.size();++i){
		if(split_FORMAT[i]=="GT"){
            GT_col=i;
            continue;
		}
        if(split_FORMAT[i]=="AD"){
            AD_col=i;
            continue;
		}
        if(split_FORMAT[i]=="DP"){
            DP_col=i;
            continue;
		}
        if(split_FORMAT[i]=="GQ"){
            GQ_col=i;
            continue;
		}
        if(split_FORMAT[i]=="PGT"){
            PGT_col=i;
            continue;
		}
        if(split_FORMAT[i]=="PID"){
            PID_col=i;
            continue;
		}
        if(split_FORMAT[i]=="PL"){
            PL_col=i;
            continue;
		}
	}
}



void AD_extractor(std::vector<std::vector<int> >& all_cols,const std::vector<std::string>& split_line,const size_t AD_col,const size_t ALT_count,const size_t FORMAT_length,const size_t FORMAT_col){
	using namespace std;
	for(size_t i(FORMAT_col+1);i<split_line.size();++i){
		const size_t II(i-FORMAT_col-1);
		const string& value(split_line[i]);
		const vector<string> split_value(split(value,':'));
		if(split_value.size()!=FORMAT_length){
			vector<int>& row(all_cols[II]);
			for(size_t j(0);j<2;++j)
				row[j]=-1;
			continue;
		}
		const string& AD_string(split_value[AD_col]);
		const vector<string> split_AD_string(split(AD_string,','));
		if(split_AD_string.size()!=ALT_count+1){
			vector<int>& row(all_cols[II]);
			for(size_t j(0);j<2;++j)
				row[j]=-1;
			continue;
		}
		
		
		vector<int>& row(all_cols[II]);
		row[0]=atoi(split_AD_string[0].c_str());
		int DPA(0);
		for(size_t j(1);j<split_AD_string.size();++j)
			DPA+=atoi(split_AD_string[j].c_str());
		row[1]=DPA;
		
	}
}






class vcf_line_cols{
	public:
	size_t CHROM_col,POS_col,ID_col,REF_col,ALT_col,QUAL_col,FILTER_col,INFO_col,FORMAT_col;
	std::vector<std::string> split_line;
	std::vector<size_t> unrel_cols;
	std::vector<size_t> trio_cols;
	size_t total_candidates;
	
	int CSQ_ExAC_AF_col;
	
	std::vector<std::vector<double> > table_L;
	
	
	
	vcf_line_cols(const std::string& vcf_header,const std::string& trio_filename,const std::string& CSQ_line){
		using namespace std;
		split_line=split(vcf_header);
		CHROM_col=POS_col=ID_col=REF_col=ALT_col=QUAL_col=FILTER_col=INFO_col=FORMAT_col=1000;
		for(size_t i(0);i<split_line.size();++i){
			std::string F(split_line[i]);
			if(F=="#CHROM")
				CHROM_col=i;
			if(F=="POS")
				POS_col=i;
			if(F=="ID")
				ID_col=i;
			if(F=="REF")
				REF_col=i;
			if(F=="ALT")
				ALT_col=i;
			if(F=="QUAL")
				QUAL_col=i;
			if(F=="FILTER")
				FILTER_col=i;
			if(F=="INFO")
				INFO_col=i;
			if(F=="FORMAT")
				FORMAT_col=i;
		}
		total_candidates=split_line.size()-FORMAT_col-1;
		
		ID_splitter(trio_filename);
		
		CSQ_ExAC_AF_col=CSQ_ExAC_AF_col_extractor(CSQ_line);
		
		
		table_L_gen(1e-5);
	}
	
	void ID_splitter(const std::string& trio_filename){
		using namespace std;
		ifstream fin(trio_filename.c_str());
		string line;
		getline(fin,line);
		const vector<string> split_line_temp(split(line));
		const string P1(split_line_temp[0]);
		const string P2(split_line_temp[1]);
		const string proband_ID(split_line_temp[2]);
		unrel_cols=vector<size_t>(total_candidates-3);
		trio_cols=vector<size_t>(3);
		size_t count(0);
		for(size_t i(0);i<total_candidates;++i){
			const size_t I(i+FORMAT_col+1);
			const string temp(split_line[I]);
			if( temp!=P1 && temp!=P2 && temp!=proband_ID ){
				unrel_cols[count]=i;
				++count;
			}
			if( temp==P1 )
				trio_cols[0]=i;
			if( temp==P2 )
				trio_cols[1]=i;
			if( temp==proband_ID )
				trio_cols[2]=i;
		}
	}
	
	int CSQ_ExAC_AF_col_extractor(const std::string& CSQ_line);
	static bool CSQ_line_checker(const std::string& line);
	
	void table_L_gen(const double mut_rate);
	
	void row_gen(std::vector<double>& row,const std::vector<size_t>& GT1,const std::vector<size_t>& GT2,const size_t& ALT_count,const double& mut_rate);
};

void vcf_line_cols::table_L_gen(const double mut_rate){
	using namespace std;
	const size_t N(1);
	const size_t combos(ALT_count_to_combos(N));
	table_L=vector<vector<double> >(combos*combos,vector<double>(combos,0));
	if(table_L.size()!=combos*combos || table_L[0].size()!=combos)
		throw(std::domain_error("error in table_gen"));
		
	for(size_t a1(0);a1<N+1;++a1){
		for(size_t a2(a1);a2<N+1;++a2){
			for(size_t a3(0);a3<N+1;++a3){
				for(size_t a4(a3);a4<N+1;++a4){
					std::vector<size_t> GT1(2,0);GT1[0]=a1;GT1[1]=a2;
					std::vector<size_t> GT2(2,0);GT2[0]=a3;GT2[1]=a4;
					const size_t I1(  N*GT1[0]+GT1[1] - int(GT1[0]*GT1[0]-GT1[0])/2  );
					const size_t I2(  N*GT2[0]+GT2[1] - int(GT2[0]*GT2[0]-GT2[0])/2  );
					const size_t II(I1*combos+I2);
					
					row_gen(table_L[II],GT1,GT2,N,mut_rate);
				}
			}
		}
	}
}

void vcf_line_cols::row_gen(std::vector<double>& row,const std::vector<size_t>& GT1,const std::vector<size_t>& GT2,const size_t& ALT_count,const double& mut_rate){
	const size_t N(ALT_count);
	const size_t combos(ALT_count_to_combos(N));
	if(row.size()!=combos)
		throw(std::domain_error("error in row_gen"));
	
	for(size_t i(0);i<combos;++i)
		row[i]=0;
	
	size_t count(0);
	for(size_t a1(0);a1<N+1;++a1){
		for(size_t a2(0);a2<N+1;++a2){
			for(size_t a3(0);a3<N+1;++a3){
				for(size_t a4(0);a4<N+1;++a4){
					double P(1);
					if(a1==GT1[0])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					if(a2==GT1[1])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					if(a3==GT2[0])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					if(a4==GT2[1])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					
					std::vector<size_t> B1(2,0);B1[0]=a1;B1[1]=a2;
					std::vector<size_t> B2(2,0);B2[0]=a3;B2[1]=a4;
					++count;
					
					for(size_t b1(0);b1<2;++b1){
						for(size_t b2(0);b2<2;++b2){
							std::vector<size_t> gt_work(2,0);
							gt_work[0]=B1[b1];
							gt_work[1]=B2[b2];
							sort(gt_work.begin(),gt_work.end());
							const size_t index(  N*gt_work[0]+gt_work[1] - int(gt_work[0]*gt_work[0]-gt_work[0])/2  );
							row[index]=row[index]+0.25*P;
						}
					}
				}
			}
		}
	}
	
	const double summ(sum(row));
	for(size_t i(0);i<row.size();++i){
		row[i]=row[i]/summ;
		row[i]=log(row[i]);
	}
}



bool vcf_line_cols::CSQ_line_checker(const std::string& line){
	using namespace std;
	if(line.size()>=14){
		const string temp(line.begin()+11,line.begin()+14);
		if(temp=="CSQ")
			return 1;
		else
			return 0;
	}
	else
		return 0;
}


int vcf_line_cols::CSQ_ExAC_AF_col_extractor(const std::string& CSQ_line){
	using namespace std;
	string::const_iterator I1,I2;
	I1=CSQ_line.begin();
	while(*I1!='"')
		++I1;
	I2=I1;
	++I2;
	while(*I2!='"')
		++I2;
	string temp(I1+1,I2);
	vector<string> split_temp( split(temp) );
	const size_t L( split_temp.size() );
	temp=split_temp[L-1];
	split_temp=split(temp,'|');
	
	for(size_t I(0);I<split_temp.size();++I){
		if(split_temp[I]=="ExAC_AF")
			return I;
	}
	return -1;
}


std::ostream& operator<<(std::ostream& out,const vcf_line_cols& obj){
	out<<"vcf_line_cols:\n**************\n";
	out<<obj.CHROM_col<<' '<<obj.POS_col<<' '<<obj.ID_col<<' '<<obj.REF_col<<' '<<obj.ALT_col<<' '<<obj.QUAL_col<<' '<<obj.FILTER_col<<' '<<obj.INFO_col<<' '<<obj.FORMAT_col<<'\n';
	out<<obj.split_line[obj.CHROM_col]<<' '<<obj.split_line[obj.POS_col]<<' '<<obj.split_line[obj.ID_col]<<' '<<obj.split_line[obj.REF_col]<<' '<<obj.split_line[obj.ALT_col]<<' '<<obj.split_line[obj.QUAL_col]<<' '<<obj.split_line[obj.FILTER_col]<<' '<<obj.split_line[obj.INFO_col]<<' '<<obj.split_line[obj.FORMAT_col]<<'\n';
	out<<"number of unrels = "<<obj.unrel_cols.size()<<"\n";
	out<<"total_candidates="<<obj.total_candidates<<"\n";
	out<<"unrel_cols=\n"<<obj.unrel_cols<<'\n';
	out<<"trio_cols:\n"<<obj.trio_cols<<"\n";
	out<<"CSQ_ExAC_AF_col="<<obj.CSQ_ExAC_AF_col<<"\n";
	out<<"table_L:\n"<<obj.table_L;
	out<<"EXP_2D(table_L):\n"<<EXP_2D(obj.table_L);
	
	return out;
}


const char *vinit[]={"intergenic_variant","feature_truncation","feature_elongation","regulatory_region_variant","regulatory_region_amplification","regulatory_region_ablation","TF_binding_site_variant","TFBS_amplification","TFBS_ablation","downstream_gene_variant","upstream_gene_variant","non_coding_transcript_variant","NMD_transcript_variant","intron_variant","non_coding_transcript_exon_variant","3_prime_UTR_variant","5_prime_UTR_variant","mature_miRNA_variant","coding_sequence_variant","synonymous_variant","stop_retained_variant","incomplete_terminal_codon_variant","splice_region_variant","missense_variant","inframe_deletion","inframe_insertion","initiator_codon_variant","stop_lost","frameshift_variant","stop_gained","splice_donor_variant","splice_acceptor_variant","transcript_ablation"};

class CSQ_data{
	public:
	std::map<std::string,double> vals;
	double ExAC_AF;
	std::string MDQ;
	int MDQ_rank;
	std::string gene;
	
	static const std::vector<std::string> ordering;
	
	CSQ_data() {}
	CSQ_data(const std::string& CSQ_line,const int ExAC_AF_col);
	
	void MDQ_update(const std::string& CON,const std::string& GENE_val);
};

const std::vector<std::string> CSQ_data::ordering(vinit,vinit+33);



void CSQ_data::MDQ_update(const std::string& CON,const std::string& GENE_val){
	using namespace std;
	vector<string> split_CON(split(CON,'&'));
	for(size_t i(0);i<split_CON.size();++i){
		const string& conseq(split_CON[i]);
		int val1(0),val2(0);
		for(int j(0);j<ordering.size();++j){
			if(MDQ==ordering[j])
				val1=j;
			if(conseq==ordering[j])
				val2=j;
		}
		if(val2>val1){
			MDQ=conseq;
			MDQ_rank=val2;
			if(GENE_val.size()>0)
				gene=GENE_val;
		}
	}
}


CSQ_data::CSQ_data(const std::string& CSQ_line,const int ExAC_AF_col): MDQ("intergenic_variant"),MDQ_rank(0) {
	using namespace std;
	const vector<string> split_line(split(CSQ_line,','));
	for(size_t i(0);i<split_line.size();++i){
		const string& each_val(split_line[i]);
		const vector<string> split_each_val(split(each_val,'|'));
		
		string freq_val;
		if(split_each_val.size()>ExAC_AF_col && ExAC_AF_col>0)
			freq_val=split_each_val[ExAC_AF_col];
		
		const string Consequence_val(split_each_val[1]);
		MDQ_update(Consequence_val,split_each_val[3]);
		if(freq_val.size()>0)
			vals[split_each_val[0]]=atof(freq_val.c_str());
		else
			vals[split_each_val[0]]=-1;
	}
	
	map<string,double>::const_iterator II(vals.begin());
	bool mark(1);
	while(II!=vals.end()){
		if(II->second>=0){
			mark=0;
			break;
		}
		++II;
	}
	if(mark){
		ExAC_AF=-1;
	}
	else{
		ExAC_AF=0;
		II=vals.begin();
		while(II!=vals.end()){
			if(II->second>=0)
				ExAC_AF+=II->second;
			++II;
		}
	}
	if(ExAC_AF>1)
		ExAC_AF=1;
}


std::ostream& operator<<(std::ostream& out,const CSQ_data& O){
	using namespace std;
	out<<"CSQ_data:"<<"\n";
	out<<"*********"<<"\n";
	map<string,double>::const_iterator O_I(O.vals.begin());
	while(O_I!=O.vals.end()){
		out<<O_I->first<<"\t"<<O_I->second<<"\n";
		++O_I;
	}
	out<<"CSQ_data::ExAC_AF="<<O.ExAC_AF<<"\n";
	out<<"CSQ_data::MDQ="<<O.MDQ<<"\n";
	out<<"CSQ_data::MDQ_rank="<<O.MDQ_rank<<"\n";
	out<<"CSQ_data::gene="<<O.gene<<"\n";
	return out;
}


class INFO_data{
	public:
	std::vector<std::string> split_INFO;
	double MQ,FS,DP,AF,ExAC_AF;
	CSQ_data CSQ_data_obj;
	
	INFO_data(){}
	INFO_data(const std::string& INFO,const int ExAC_AF_col);
};

INFO_data::INFO_data(const std::string& INFO,const int ExAC_AF_col): split_INFO(split(INFO,';')) ,MQ(-1),FS(-1),DP(-1),AF(-1),ExAC_AF(-1) {
	std::vector<std::string>::const_iterator i(split_INFO.begin());
	for(;i!=split_INFO.end();++i){
		const std::string& temp(*i);
		if(temp[0]=='M' && temp[1]=='Q' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			MQ=atof(new_temp.c_str());
			continue;
		}
		if(temp[0]=='F' && temp[1]=='S' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			FS=atof(new_temp.c_str());
			continue;
		}
		if(temp[0]=='D' && temp[1]=='P' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			DP=atof(new_temp.c_str());
			continue;
		}
		if(temp[0]=='A' && temp[1]=='F' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			AF=atof(new_temp.c_str());
			continue;
		}
		if(std::string(temp.begin(),temp.begin()+8)=="ExAC_AF="){
			const std::string new_temp(temp.begin()+8,temp.end());
			ExAC_AF=atof(new_temp.c_str());
		}
		if(std::string(temp.begin(),temp.begin()+4)=="CSQ="){
			const std::string new_temp(temp.begin()+4,temp.end());
			CSQ_data_obj=CSQ_data(new_temp,ExAC_AF_col);
		}
	}
}
std::ostream& operator<<(std::ostream& out,const INFO_data& O){
	out<<"INFO_data:\n**********\n";
	out<<"MQ="<<O.MQ<<'\n';
	out<<"FS="<<O.FS<<'\n';
	out<<"DP="<<O.DP<<'\n';
	out<<"AF="<<O.AF<<'\n';
	out<<"ExAC_AF="<<O.ExAC_AF<<"\n";
	out<<O.CSQ_data_obj;
	return out;
}


std::ostream& operator<<(std::ostream& out,const INFO_data& O);

class location{
	public:
	size_t chrom,pos;
	
	location(const size_t& C,const size_t& P): chrom(C),pos(P) {}
	location(){}
	
};
bool operator<(const location& loc1,const location& loc2);
bool operator==(const location& loc1,const location& loc2);
bool operator>(const location& loc1,const location& loc2);
bool operator<=(const location& loc1,const location& loc2);
bool operator>=(const location& loc1,const location& loc2);
bool operator!=(const location& loc1,const location& loc2);
long operator-(const location& L1,const location& L2);


std::ostream& operator<<(std::ostream& out,const location& loc){
	out<<loc.chrom<<"\t"<<loc.pos;
	return out;
}

size_t str_to_size_t(const std::string& S);


void AD_extractor_2(std::vector<std::vector<int> >& AD_list,std::istringstream& is,const size_t AD_col,const size_t ALT_count,const size_t FORMAT_length,const size_t FORMAT_col){
	using namespace std;
	const int LL(AD_list.size());
	string val;
	for(int i(0);i<LL;++i){
		is>>val;
		vector<string> split_val( split(val,':') );
		vector<int>& AD(AD_list[i]);
		if(split_val.size() != FORMAT_length){
			AD[0]=-1;
			AD[1]=-1;
			continue;
		}
		
		const string& AD_string(split_val[AD_col]);
		vector<string> split_AD_string( split(AD_string,',') );
		if(split_AD_string.size() != ALT_count+1){
			AD[0]=-1;
			AD[1]=-1;
			continue;
		}
		
		AD[0]=atoi(split_AD_string[0].c_str());
		AD[1]=0;
		for(int j(1);j<split_AD_string.size();++j){
			AD[1] += atoi(split_AD_string[j].c_str());
		}
	}
}




class each_line_data_2{
	// an object of this class is declred once and mem_set is used after that
	public:
	std::vector<std::string> split_line;
	location LOC;
	std::string REF;
	std::string ALT;
	std::vector<std::string> split_ALT;
	const size_t combos;
	size_t ALT_count,combos_actual;
	std::string FORMAT;
	std::vector<std::string> split_FORMAT;
	size_t GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col;
	std::string INFO;
	size_t FORMAT_length;
	std::vector<std::vector<int> > AD_list;
	std::vector<std::vector<double> > GL_L_list;
	INFO_data INFO_data_obj;
	
	
	
	
	each_line_data_2(const std::string& line,const vcf_line_cols& obj);
	
	void mem_set  (const std::string& line,const vcf_line_cols& obj);
	
	double HET_P_calc(const vcf_line_cols& VLCO) const;
};

double each_line_data_2::HET_P_calc(const vcf_line_cols& VLCO) const{
	using namespace std;
	
	const size_t person_I( VLCO.trio_cols[2] );
	const size_t P1_I( VLCO.trio_cols[0] );
	const size_t P2_I( VLCO.trio_cols[1] );
	
	const vector<int>& AD( AD_list[person_I] );
	const vector<int>& AD_P1( AD_list[P1_I] );
	const vector<int>& AD_P2( AD_list[P2_I] );
	
	const vector<double>& GL_L( GL_L_list[person_I] );
	const vector<double>& GL_L_P1( GL_L_list[P1_I] );
	const vector<double>& GL_L_P2( GL_L_list[P2_I] );
	
	if(AD[0]>=0 && AD[1]>=0 && AD_P1[0]>=0 && AD_P1[1]>=0 && AD_P2[0]>=0 && AD_P2[1]>=0){
		vector<vector<double> > work_table( VLCO.table_L );
		double max_val( -numeric_limits<double>::infinity() );
		for(size_t I1(0);I1<3;++I1){
			for(size_t I2(0);I2<3;++I2){
				for(size_t I3(0);I3<3;++I3){
					const size_t index(3*I1+I2);
					work_table[index][I3]+=( GL_L_P1[I1]+GL_L_P2[I2]+GL_L[I3] );
					if( work_table[index][I3]>max_val )
						max_val=work_table[index][I3];
				}
			}
		}
		double summ(0);
		double summ2(0);
		for(size_t I1(0);I1<3;++I1){
			for(size_t I2(0);I2<3;++I2){
				for(size_t I3(0);I3<3;++I3){
					const size_t index(3*I1+I2);
					work_table[index][I3]-=max_val;
					work_table[index][I3]=exp( work_table[index][I3] );
					summ+=work_table[index][I3];
					if( I3==1 && index>=1 && index<=7 || I3==2 && ( index>=4 && index<=5 || index>=7 && index<=8 ) )
					//if( I3==1 && ( index==1 || index==3 || index==4 ) )
						summ2+=work_table[index][I3];
				}
			}
		}
		
		//work_table=work_table*(1/summ);
		//cout<<"work_table:\n"<<work_table;
		
		return ( summ2/summ );
	}
	else
		return 0;
}



each_line_data_2::each_line_data_2(const std::string& line,const vcf_line_cols& obj): combos(3) {
	using namespace std;
	split_line=split(line);
	const size_t C(str_to_size_t(split_line[obj.CHROM_col]));
	const size_t P(atoi(split_line[obj.POS_col].c_str()));
	LOC=location(C,P);
	REF=split_line[obj.REF_col];
	ALT=split_line[obj.ALT_col];
	split_ALT=split(ALT,',');
	ALT_count=split_ALT.size();
	combos_actual=ALT_count_to_combos(ALT_count);
	FORMAT=split_line[obj.FORMAT_col];
	split_FORMAT=split(FORMAT,':');
	col_extractor_for_FORMAT(FORMAT,GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col);
	INFO=split_line[obj.INFO_col];
	FORMAT_length=split_FORMAT.size();
	AD_list  =vector<vector<int> >(obj.total_candidates,vector<int>(2,-1));
	AD_extractor(AD_list,split_line,AD_col,ALT_count,FORMAT_length,obj.FORMAT_col);
	GL_L_list=vector<vector<double> >(obj.total_candidates,vector<double>(3,100));
	INFO_data_obj=INFO_data(INFO,obj.CSQ_ExAC_AF_col);
}

void each_line_data_2::mem_set(const std::string& line,const vcf_line_cols& obj){
	using namespace std;
	istringstream is(line);
	string temp;
	is>>temp;const size_t C(str_to_size_t(temp));
	is>>temp;const size_t P(atoi(temp.c_str()));
	LOC=location(C,P);
	is>>temp;
	is>>REF>>ALT;
	split_ALT=split(ALT,',');
	ALT_count=split_ALT.size();
	combos_actual=ALT_count_to_combos(ALT_count);
	is>>temp;
	is>>temp;
	is>>INFO;
	is>>FORMAT;
	split_FORMAT=split(FORMAT,':');
	col_extractor_for_FORMAT(FORMAT,GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col);
	FORMAT_length=split_FORMAT.size();
	//AD_list  =vector<vector<int> >(obj.total_candidates,vector<int>(2,-1));
	AD_extractor_2(AD_list,is,AD_col,ALT_count,FORMAT_length,obj.FORMAT_col);
	//GL_L_list=vector<vector<double> >(obj.total_candidates,vector<double>(3,100));
	INFO_data_obj=INFO_data(INFO,obj.CSQ_ExAC_AF_col);
}

std::ostream& operator<<(std::ostream& out,const each_line_data_2& obj){
	out<<"each_line_data_2:\n***************\n";
	out<<"LOC="<<obj.LOC<<"\n";
	out<<"REF="<<obj.REF<<"\n";
	out<<"ALT="<<obj.ALT<<'\n';
	out<<"split_ALT="<<obj.split_ALT<<'\n';
	out<<"ALT_count, combos,combos_actual = "<<obj.ALT_count<<", "<<obj.combos<<", "<<obj.combos_actual<<'\n';
	out<<"GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col="<<obj.GT_col<<' '<<obj.AD_col<<' '<<obj.DP_col<<' '<<obj.GQ_col<<' '<<obj.PGT_col<<' '<<obj.PID_col<<' '<<obj.PL_col<<'\n';
	out<<"FORMAT="<<obj.FORMAT<<'\n';
	out<<"FORMAT_length="<<obj.FORMAT_length<<'\n';
	out<<"INFO=  "<<obj.INFO<<'\n';
	out<<"INFO_data_obj=\n"<<obj.INFO_data_obj;
	out<<"AD_list=\n";
	for(size_t i(0);i<obj.AD_list.size();++i)
		out<<i<<")\t"<<obj.AD_list[i]<<"\t"<<obj.GL_L_list[i]<<'\n';
	out<<"AD_list.size()="<<obj.AD_list.size()<<'\n';
	//out<<"max_sum_data_obj=\n"<<obj.max_sum_data_obj;
	
	return out;
}

class each_line_data{
	// an object of this class is declred once and mem_set is used after that
	public:
	std::vector<std::string> split_line;
	location LOC;
	std::string REF;
	std::string ALT;
	std::vector<std::string> split_ALT;
	const size_t combos;
	size_t ALT_count,combos_actual;
	std::string FORMAT;
	std::vector<std::string> split_FORMAT;
	size_t GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col;
	std::string INFO;
	size_t FORMAT_length;
	std::vector<std::vector<int> > AD_list;
	std::vector<std::vector<double> > GL_L_list;
	INFO_data INFO_data_obj;
	
	
	
	
	each_line_data(const std::string& line,const vcf_line_cols& obj);
	
	void mem_set  (const std::string& line,const vcf_line_cols& obj);
	
	double HET_P_calc(const vcf_line_cols& VLCO) const;
};

double each_line_data::HET_P_calc(const vcf_line_cols& VLCO) const{
	using namespace std;
	
	const size_t person_I( VLCO.trio_cols[2] );
	const size_t P1_I( VLCO.trio_cols[0] );
	const size_t P2_I( VLCO.trio_cols[1] );
	
	const vector<int>& AD( AD_list[person_I] );
	const vector<int>& AD_P1( AD_list[P1_I] );
	const vector<int>& AD_P2( AD_list[P2_I] );
	
	const vector<double>& GL_L( GL_L_list[person_I] );
	const vector<double>& GL_L_P1( GL_L_list[P1_I] );
	const vector<double>& GL_L_P2( GL_L_list[P2_I] );
	
	if(AD[0]>=0 && AD[1]>=0 && AD_P1[0]>=0 && AD_P1[1]>=0 && AD_P2[0]>=0 && AD_P2[1]>=0){
		vector<vector<double> > work_table( VLCO.table_L );
		double max_val( -numeric_limits<double>::infinity() );
		for(size_t I1(0);I1<3;++I1){
			for(size_t I2(0);I2<3;++I2){
				for(size_t I3(0);I3<3;++I3){
					const size_t index(3*I1+I2);
					work_table[index][I3]+=( GL_L_P1[I1]+GL_L_P2[I2]+GL_L[I3] );
					if( work_table[index][I3]>max_val )
						max_val=work_table[index][I3];
				}
			}
		}
		double summ(0);
		double summ2(0);
		for(size_t I1(0);I1<3;++I1){
			for(size_t I2(0);I2<3;++I2){
				for(size_t I3(0);I3<3;++I3){
					const size_t index(3*I1+I2);
					work_table[index][I3]-=max_val;
					work_table[index][I3]=exp( work_table[index][I3] );
					summ+=work_table[index][I3];
					if( I3==1 && index>=1 && index<=7 || I3==2 && ( index>=4 && index<=5 || index>=7 && index<=8 ) )
					//if( I3==1 && ( index==1 || index==3 || index==4 ) )
						summ2+=work_table[index][I3];
				}
			}
		}
		
		//work_table=work_table*(1/summ);
		//cout<<"work_table:\n"<<work_table;
		
		return ( summ2/summ );
	}
	else
		return 0;
}


each_line_data::each_line_data(const std::string& line,const vcf_line_cols& obj): combos(3) {
	using namespace std;
	split_line=split(line);
	const size_t C(str_to_size_t(split_line[obj.CHROM_col]));
	const size_t P(atoi(split_line[obj.POS_col].c_str()));
	LOC=location(C,P);
	REF=split_line[obj.REF_col];
	ALT=split_line[obj.ALT_col];
	split_ALT=split(ALT,',');
	ALT_count=split_ALT.size();
	combos_actual=ALT_count_to_combos(ALT_count);
	FORMAT=split_line[obj.FORMAT_col];
	split_FORMAT=split(FORMAT,':');
	col_extractor_for_FORMAT(FORMAT,GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col);
	INFO=split_line[obj.INFO_col];
	FORMAT_length=split_FORMAT.size();
	AD_list  =vector<vector<int> >(obj.total_candidates,vector<int>(2,-1));
	AD_extractor(AD_list,split_line,AD_col,ALT_count,FORMAT_length,obj.FORMAT_col);
	GL_L_list=vector<vector<double> >(obj.total_candidates,vector<double>(3,100));
	INFO_data_obj=INFO_data(INFO,obj.CSQ_ExAC_AF_col);
}

void each_line_data::mem_set(const std::string& line,const vcf_line_cols& obj){
	using namespace std;
	split_line=split(line);
	const size_t C(str_to_size_t(split_line[obj.CHROM_col]));
	const size_t P(atoi(split_line[obj.POS_col].c_str()));
	LOC=location(C,P);
	REF=split_line[obj.REF_col];
	ALT=split_line[obj.ALT_col];
	split_ALT=split(ALT,',');
	ALT_count=split_ALT.size();
	combos_actual=ALT_count_to_combos(ALT_count);
	FORMAT=split_line[obj.FORMAT_col];
	split_FORMAT=split(FORMAT,':');
	col_extractor_for_FORMAT(FORMAT,GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col);
	INFO=split_line[obj.INFO_col];
	FORMAT_length=split_FORMAT.size();
	//AD_list  =vector<vector<int> >(obj.total_candidates,vector<int>(2,-1));
	AD_extractor(AD_list,split_line,AD_col,ALT_count,FORMAT_length,obj.FORMAT_col);
	//GL_L_list=vector<vector<double> >(obj.total_candidates,vector<double>(3,100));
	INFO_data_obj=INFO_data(INFO,obj.CSQ_ExAC_AF_col);
}



std::ostream& operator<<(std::ostream& out,const each_line_data& obj){
	out<<"each_line_data:\n***************\n";
	out<<"LOC="<<obj.LOC<<"\n";
	out<<"REF="<<obj.REF<<"\n";
	out<<"ALT="<<obj.ALT<<'\n';
	out<<"split_ALT="<<obj.split_ALT<<'\n';
	out<<"ALT_count, combos,combos_actual = "<<obj.ALT_count<<", "<<obj.combos<<", "<<obj.combos_actual<<'\n';
	out<<"GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col="<<obj.GT_col<<' '<<obj.AD_col<<' '<<obj.DP_col<<' '<<obj.GQ_col<<' '<<obj.PGT_col<<' '<<obj.PID_col<<' '<<obj.PL_col<<'\n';
	out<<"FORMAT="<<obj.FORMAT<<'\n';
	out<<"FORMAT_length="<<obj.FORMAT_length<<'\n';
	out<<"INFO=  "<<obj.INFO<<'\n';
	out<<"INFO_data_obj=\n"<<obj.INFO_data_obj;
	out<<"AD_list=\n";
	for(size_t i(0);i<obj.AD_list.size();++i)
		out<<i<<")\t"<<obj.AD_list[i]<<"\t"<<obj.GL_L_list[i]<<'\n';
	out<<"AD_list.size()="<<obj.AD_list.size()<<'\n';
	//out<<"max_sum_data_obj=\n"<<obj.max_sum_data_obj;
	
	return out;
}









bool operator<(const location& loc1,const location& loc2){
	if(loc1.chrom<loc2.chrom)
		return 1;
	else{
		if(loc1.chrom==loc2.chrom){
			if(loc1.pos<loc2.pos)
				return 1;
			else
				return 0;
		}
		else
			return 0;
	}
}

bool operator==(const location& loc1,const location& loc2){
	if(loc1.chrom==loc2.chrom && loc1.pos==loc2.pos)
		return 1;
	else
		return 0;
}

bool operator!=(const location& loc1,const location& loc2){
	return !(loc1==loc2);
}

bool operator>(const location& loc1,const location& loc2){
	if( !(loc1<loc2) && !(loc1==loc2) )
		return 1;
	else
		return 0;
}

bool operator<=(const location& loc1,const location& loc2){
	if(loc1<loc2 || loc1==loc2)
		return 1;
	else
		return 0;
}

bool operator>=(const location& loc1,const location& loc2){
	if(loc1>loc2 || loc1==loc2)
		return 1;
	else
		return 0;
}

long operator-(const location& L1,const location& L2){
	if(L1.chrom==L2.chrom)
		return long(L1.pos)-long(L2.pos);
	else{
		if(L1.chrom>L2.chrom)
			return std::numeric_limits<long>::max();
		else
			return -std::numeric_limits<long>::max();
	}
}





size_t str_to_size_t(const std::string& S){
	std::string S2,S_temp;
	if(S.size()>3)
		S_temp=std::string(S.begin(),S.begin()+3);
	else
		S_temp=S;
	
	if(S_temp=="chr")
		S2=std::string(S.begin()+3,S.end());
	else
		S2=S;
	
	//std::cout<<"S2="<<S2<<'\n';
	
	const size_t num(atoi(S2.c_str()));
	if(num>=1 && num<=22)
		return num;
	else
		if(S2==std::string("M") || S2==std::string("m"))
			return 0;
		else
			if(S2==std::string("X") || S2==std::string("x"))
				return 23;
			else
				if(S2==std::string("Y") || S2==std::string("y"))
					return 24;
				else
					return 25;
}
///////
class parameters2{
	public:
	double rho_new,rho_old;
	std::vector<double> GF_new,GF_old,GF_L_new,GF_L_old;
	double AF_new,AF_old;
	double log_likelihood_new,log_likelihood_old;
	size_t extra_iter;
	double AF_unrel;
	
	parameters2();
	void reset();
	
	void update_GL_L_ALL(const vcf_line_cols& VLCO,each_line_data& ELDO);
	
	void EM_step(const vcf_line_cols&,each_line_data&);
	
	void EM_full(const vcf_line_cols&,each_line_data&);
	
	void AF_unrel_calc(const vcf_line_cols& VLCO,each_line_data& ELDO);
};

std::ostream& operator<<(std::ostream& out,const parameters2& O){
	using namespace std;
	out<<"parameters2:\n";
	out<<"************\n";
	out<<"rho new and old:\n";
	out<<O.rho_new<<"\t"<<O.rho_old<<"\n";
	out<<"GF new and old:\n";
	out<<O.GF_new<<"\n"<<O.GF_old<<"\n";
	out<<"GF_L new and old:\n";
	out<<O.GF_L_new<<"\n"<<O.GF_L_old<<"\n";
	out<<"AF new and old:\n";
	out<<O.AF_new<<"\t"<<O.AF_old<<"\n";
	out<<"log_likelihood new and old:\n";
	out<<O.log_likelihood_new<<"\t"<<O.log_likelihood_old<<"\n";
	out<<"extra_iter="<<O.extra_iter<<"\n";
	out<<"AF_unrel="<<O.AF_unrel<<"\n";
}


parameters2::parameters2(){
	using namespace std;
	rho_new=0.8;
	rho_old=0.8;
	AF_new=0.1;
	AF_old=0.1;
	GF_new=vector<double>(3);
	GF_new[0]=0.5;
	GF_new[1]=0.3;
	GF_new[2]=0.2;
	GF_old=GF_new;
	GF_L_new=LOG_1D(GF_new);
	GF_L_old=GF_L_new;
	AF_new=GF_new[1]/2.0+GF_new[2];
	AF_old=AF_new;
	
	log_likelihood_old=0;
	log_likelihood_new=0;
	extra_iter=0;
	AF_unrel=-1;
}

void parameters2::reset(){
	using namespace std;
	rho_new=0.8;
	rho_old=0.8;
	AF_new=0.1;
	AF_old=0.1;
	GF_new=vector<double>(3);
	GF_new[0]=0.5;
	GF_new[1]=0.3;
	GF_new[2]=0.2;
	GF_old=GF_new;
	GF_L_new=LOG_1D(GF_new);
	GF_L_old=GF_L_new;
	AF_new=GF_new[1]/2.0+GF_new[2];
	AF_old=AF_new;
	
	log_likelihood_old=0;
	log_likelihood_new=0;
	extra_iter=0;
	AF_unrel=-1;
}


void parameters2::EM_step(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	double ExAC_AF( ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF );
	if( ExAC_AF<0 )
		ExAC_AF=0;
	
	
	///////
	rho_old=rho_new;
	GF_old=GF_new;
	GF_L_old=GF_L_new;
	AF_old=AF_new;
	log_likelihood_old=log_likelihood_new;
	///////
	
	double T1(1),T2(1);
	vector<double> FT(3,1);
	double NN(0);
	
	const vector<size_t>& unrel_cols(VLCO.unrel_cols);
	const vector<vector<int> >& AD_list( ELDO.AD_list );
	vector<vector<double> >& GL_L_list( ELDO.GL_L_list );
	
	const double rho_L( log(rho_old) );
	const double rho_L_1( log(1-rho_old) );
	const double mid( log(0.5) );
	
	log_likelihood_new=rho_L+rho_L_1+sum(GF_L_old);
	for(size_t i(0);i<unrel_cols.size();++i){
		const size_t index( unrel_cols[i] );
		const vector<int>& AD( AD_list[ index ] );
		vector<double>& GL_L( GL_L_list[ index ] );
		if(AD[0]>=0 && AD[1]>=0 && !(AD[0]==0 && AD[1]==0) ){
			GL_L[0]=AD[0]*rho_L+AD[1]*rho_L_1;
			GL_L[1]=(AD[0]+AD[1])*mid;
			GL_L[2]=AD[1]*rho_L+AD[0]*rho_L_1;
			vector<double> PP( GL_L+GF_L_old );
			const double max_val(max(PP));
			PP=PP+(-max_val);
			PP=EXP_1D(PP);
			const double summ( sum(PP) );
			log_likelihood_new+=log(summ)+max_val;
			PP=PP*(1/summ);
			T1+=AD[0]*PP[0]+AD[1]*PP[2];
			T2+=AD[1]*PP[0]+AD[0]*PP[2];
			
			if(FT.size()!=3)
				throw(domain_error("ERROR1 in parameters2::EM_step"));
			
			FT[0]+=PP[0];
			FT[1]+=PP[1];
			FT[2]+=PP[2];
			
			++NN;
		}
	}
	
	
	rho_new=1/(1+T2/T1);
	GF_new= FT*(1/(NN+3));
	GF_L_new=LOG_1D(GF_new);
	AF_new=GF_new[1]/2+GF_new[2];
	
	
	
	//cout.precision(30);
	//cout<<log_likelihood_new<<",";
}


void parameters2::EM_full(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	reset();
	for(size_t i(0);i<3;++i)
		EM_step(VLCO,ELDO);
	extra_iter=0;
	const size_t max_iter(500);
	while( abs(log_likelihood_new-log_likelihood_old)>1e-10 ){
		EM_step(VLCO,ELDO);
		++extra_iter;
		 if(extra_iter>=max_iter)
			break;
	}
	AF_unrel_calc(VLCO,ELDO);
}



void parameters2::AF_unrel_calc(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	
	const double L_rho(log(rho_old));
	const double L_1_rho(log(1-rho_old));
	const double mid(log(0.5));
	const vector<size_t>& unrel_cols(VLCO.unrel_cols);
	const vector<vector<int> >& AD_list(ELDO.AD_list);
	vector<double> GL_L(3,0);
	AF_unrel=0;
	double denominator(0);
	for(size_t i(0);i<unrel_cols.size();++i){
		const size_t I(unrel_cols[i]);
		const vector<int>& AD(AD_list[I]);
		if(AD[0]>=0 && AD[1]>=0){
			GL_L[0]=AD[0]*L_rho+AD[1]*L_1_rho  + GF_L_old[0];
			GL_L[1]=(AD[0]+AD[1])*mid          + GF_L_old[1];
			GL_L[2]=AD[1]*L_rho+AD[0]*L_1_rho  + GF_L_old[2];
			const double max_val(max(GL_L));
			double summ(0);
			for(size_t j(0);j<GL_L.size();++j){
				GL_L[j]=GL_L[j]-max_val;
				GL_L[j]=exp(GL_L[j]);
				summ+=GL_L[j];
			}
			for(size_t j(0);j<GL_L.size();++j){
				GL_L[j]=GL_L[j]/summ;
			}
			AF_unrel+= ( 0.5*GL_L[1]+GL_L[2] );
			denominator+=1;
		}
	}
	AF_unrel=AF_unrel/denominator;
}




void parameters2::update_GL_L_ALL(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	const double L_rho(log(rho_old));
	const double L_1_rho(log(1-rho_old));
	const double mid(log(0.5));
	//const vector<size_t>& unrel_cols(VLCO.unrel_cols);
	const vector<vector<int> >& AD_list(ELDO.AD_list);
	vector<vector<double> >& GL_L_list(ELDO.GL_L_list);
	if(AD_list.size()!=GL_L_list.size())
		throw(domain_error("ERROR in parameters2::update_GL_L_ALL"));
	
	for(size_t i(0);i<AD_list.size();++i){
		const vector<int>& AD(AD_list[i]);
		vector<double>& GL_L(GL_L_list[i]);
		if(AD[0]>=0 && AD[1]>=0){
			GL_L[0]=AD[0]*L_rho+AD[1]*L_1_rho;
			GL_L[1]=(AD[0]+AD[1])*mid;
			GL_L[2]=AD[1]*L_rho+AD[0]*L_1_rho;
			/*
			const double max_val(max(GL_L));
			for(size_t j(0);j<GL_L.size();++j){
				GL_L[j]=GL_L[j]-max_val;
			}
			*/
		}
	}
}





///////

/*
class parameters{
	public:
	double rho_new,rho_old,AF_new,AF_old;
	std::vector<double> GF_L_new,GF_L_old;
	double log_likelihood_old,log_likelihood_new;
	size_t extra_iter;
	double AF_unrel;
	
	parameters();
	
	void update_GL_L_ALL(const vcf_line_cols& VLCO,each_line_data& ELDO);
	
	void EM_step(const vcf_line_cols& VLCO,each_line_data& ELDO);
	
	void EM_full(const vcf_line_cols& VLCO,each_line_data& ELDO);
	
	void AF_unrel_calc(const vcf_line_cols& VLCO,each_line_data& ELDO);
};


void parameters::AF_unrel_calc(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	const double GFF_L[]={    2*log(1-AF_old)  ,  log(2)+log(AF_old)+log(1-AF_old)  ,  2*log(AF_old)    };
	const vector<double> GF_L(GFF_L,GFF_L+3);
	if(GF_L.size()!=3)
		throw(domain_error("ERROR1 in parameters::AF_unrel_calc"));
	
	const double L_rho(log(rho_old));
	const double L_1_rho(log(1-rho_old));
	const double mid(log(0.5));
	const vector<size_t>& unrel_cols(VLCO.unrel_cols);
	const vector<vector<int> >& AD_list(ELDO.AD_list);
	vector<vector<double> >& GL_L_list(ELDO.GL_L_list);
	AF_unrel=0;
	double denominator(0);
	for(size_t i(0);i<unrel_cols.size();++i){
		const size_t I(unrel_cols[i]);
		const vector<int>& AD(AD_list[I]);
		vector<double>& GL_L(GL_L_list[I]);
		if(AD[0]>=0 && AD[1]>=0){
			GL_L[0]=AD[0]*L_rho+AD[1]*L_1_rho  + GF_L[0];
			GL_L[1]=(AD[0]+AD[1])*mid          + GF_L[1];
			GL_L[2]=AD[1]*L_rho+AD[0]*L_1_rho  + GF_L[2];
			const double max_val(max(GL_L));
			double summ(0);
			for(size_t j(0);j<GL_L.size();++j){
				GL_L[j]=GL_L[j]-max_val;
				GL_L[j]=exp(GL_L[j]);
				summ+=GL_L[j];
			}
			for(size_t j(0);j<GL_L.size();++j){
				GL_L[j]=GL_L[j]/summ;
			}
			AF_unrel+= ( 0.5*GL_L[1]+GL_L[2] );
			denominator+=1;
		}
	}
	AF_unrel=AF_unrel/denominator;
}



parameters::parameters(){
	using namespace std;
	rho_new=0.8;
	rho_old=0.8;
	AF_new=0.1;
	AF_old=0.1;
	GF_L_new=vector<double>(3);
	
	GF_L_new[0]=2*log(1-AF_new);
	GF_L_new[1]=log(2)+log(AF_new)+log(1-AF_new);
	GF_L_new[2]=2*log(AF_new);
	//GF_L_new=GF_L_new+(-max(GF_L_new));
	
	GF_L_old=GF_L_new;
	
	log_likelihood_old=0;
	log_likelihood_new=0;
	extra_iter=0;
	AF_unrel=-1;
}


void parameters::update_GL_L_ALL(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	const double L_rho(log(rho_old));
	const double L_1_rho(log(1-rho_old));
	const double mid(log(0.5));
	const vector<size_t>& unrel_cols(VLCO.unrel_cols);
	const vector<vector<int> >& AD_list(ELDO.AD_list);
	vector<vector<double> >& GL_L_list(ELDO.GL_L_list);
	for(size_t i(0);i<AD_list.size();++i){
		const size_t I(i);
		const vector<int>& AD(AD_list[I]);
		vector<double>& GL_L(GL_L_list[I]);
		if(AD[0]>=0 && AD[1]>=0){
			GL_L[0]=AD[0]*L_rho+AD[1]*L_1_rho;
			GL_L[1]=(AD[0]+AD[1])*mid;
			GL_L[2]=AD[1]*L_rho+AD[0]*L_1_rho;
			
			const double max_val(max(GL_L));
			for(size_t j(0);j<GL_L.size();++j){
				GL_L[j]=GL_L[j]-max_val;
			}
			
		}
	}
}




void parameters::EM_full(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	
	
	//ofstream fout("to_plot.txt");
	
	for(size_t i(0);i<3;i++){
		EM_step(VLCO,ELDO);
		//fout<<log_likelihood_new<<"\t";
	}
	extra_iter=0;
	const size_t max_iter(300);
	while( (log_likelihood_new-log_likelihood_old)>1e-8 ){
		EM_step(VLCO,ELDO);
		//fout<<log_likelihood_new<<"\t";
		++extra_iter;
		
		
		
		if(extra_iter>=max_iter)
			break;
		
	}
	
	AF_unrel_calc(VLCO,ELDO);
}





void parameters::EM_step(const vcf_line_cols& VLCO,each_line_data& ELDO){
	using namespace std;
	//update_GL_L(VLCO,ELDO);
	
	double ExAC_AF( ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF );
	if( ExAC_AF<0 )
		ExAC_AF=0;
	
	//double T1_AF(ExAC_AF*1000+1.0-1.0),T2_AF( (1.0-ExAC_AF)*1000+1.0-1.0 ),T1_rho(1),T2_rho(1);
	double T1_AF(1),T2_AF(1),T1_rho(1),T2_rho(1);
	
	const double L_rho(log(rho_old));
	const double L_1_rho(log(1-rho_old));
	const double mid(log(0.5));
	
	const vector<size_t>& unrel_cols(VLCO.unrel_cols);
	const vector<vector<int> >& AD_list(ELDO.AD_list);
	vector<vector<double> >& GL_L_list(ELDO.GL_L_list);
	rho_old=rho_new;
	AF_old=AF_new;
	GF_L_old=GF_L_new;
	log_likelihood_old=log_likelihood_new;
	
	log_likelihood_new=L_rho+L_1_rho+log(AF_old)+log(1-AF_old);
	for(size_t i(0);i<unrel_cols.size();++i){
		const size_t I(unrel_cols[i]);
		const vector<int>& AD(AD_list[I]);
		vector<double>& GL_L(GL_L_list[I]);
		if(AD[0]>=0 && AD[1]>=0){
			GL_L[0]=AD[0]*L_rho+AD[1]*L_1_rho;
			GL_L[1]=(AD[0]+AD[1])*mid;
			GL_L[2]=AD[1]*L_rho+AD[0]*L_1_rho;
			vector<double> PP(GL_L+GF_L_old);
			const double max_val(max(PP));
			PP=PP+(-max_val);
			PP=EXP_1D(PP);
			const double summ(sum(PP));
			log_likelihood_new+=log(summ)+max_val;
			PP=PP*(1/summ);
			T1_AF+=PP[1]+2*PP[2];
			T2_AF+=2*PP[0]+PP[1];
			T1_rho+=AD[0]*PP[0]+AD[1]*PP[2];
			T2_rho+=AD[1]*PP[0]+AD[0]*PP[2];
		}
	}
	
	
	rho_new=1.0/(1.0+T2_rho/T1_rho);
	
	if(rho_new<0.5)
		rho_new=0.5;
	
	AF_new=1.0/(1.0+T2_AF/T1_AF);
	GF_L_new[0]=2*log(1-AF_new);
	GF_L_new[1]=log(2)+log(AF_new)+log(1-AF_new);
	GF_L_new[2]=2*log(AF_new);
	
}





std::ostream& operator<<(std::ostream& out,const parameters& O){
	out<<"parameters:\n";
	out<<"***********\n";
	out<<"AFs="<<O.AF_new<<"\t"<<O.AF_old<<"\n";
	out<<"rhos="<<O.rho_new<<"\t"<<O.rho_old<<"\n";
	out<<"GF_L_new="<<O.GF_L_new<<"\n";
	out<<"GF_L_old="<<O.GF_L_old<<"\n";
	out<<"log_likelihood_new="<<O.log_likelihood_new<<"\n";
	out<<"log_likelihood_old="<<O.log_likelihood_old<<"\n";
	out<<"extra_iter="<<O.extra_iter<<"\n";
	out<<"AF_unrel="<<O.AF_unrel<<"\n";
	return out;
}

*/

class call_data{
	public:
	double Probty;
	location LOC;
	std::vector<double> GL_L;
	std::vector<int> AD;
	std::vector<double> GL_L_P1,GL_L_P2;
	std::vector<int> AD_P1,AD_P2;
	double ExAC_AF,AF_unrel;
	int MDQ_rank;
	std::string MDQ;
	std::string gene;
	
	call_data(const vcf_line_cols& VLCO,const each_line_data& ELDO,const parameters2& P,const double probty);
	
	//void binary_write(std::ofstream& fout) const;
	//void binary_read (std::ifstream& fin);
};






call_data::call_data(const vcf_line_cols& VLCO,const each_line_data& ELDO,const parameters2& P,const double probty){
	using namespace std;
	Probty=probty;
	LOC=ELDO.LOC;
	const size_t person_I( VLCO.trio_cols[2] );
	const size_t P1_I( VLCO.trio_cols[0] );
	const size_t P2_I( VLCO.trio_cols[1] );
	const vector<vector<double> >& GL_L_list( ELDO.GL_L_list );
	const vector<vector<int> >& AD_list( ELDO.AD_list );
	
	AD=AD_list[person_I];
	AD_P1=AD_list[P1_I];
	AD_P2=AD_list[P2_I];
	
	GL_L=GL_L_list[person_I];
	GL_L_P1=GL_L_list[P1_I];
	GL_L_P2=GL_L_list[P2_I];
	
	ExAC_AF=ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF;
	AF_unrel=P.AF_unrel;
	MDQ_rank=ELDO.INFO_data_obj.CSQ_data_obj.MDQ_rank;
	MDQ=ELDO.INFO_data_obj.CSQ_data_obj.MDQ;
	gene=ELDO.INFO_data_obj.CSQ_data_obj.gene;
}


std::ostream& operator<<(std::ostream& out,const call_data& O){
	using namespace std;
	out<<"call_data:\n";
	out<<"**********\n";
	out<<"Probty="<<O.Probty<<"\n";
	out<<"LOC="<<O.LOC<<"\n";
	out<<"ADs:\n";
	out<<O.AD<<"\n"<<O.AD_P1<<"\n"<<O.AD_P2<<"\n";
	out<<"GL_Ls:\n";
	out<<O.GL_L<<"\n"<<O.GL_L_P1<<"\n"<<O.GL_L_P2<<"\n";
	out<<"ExAC_AF="<<O.ExAC_AF<<"\n";
	out<<"AF_unrel="<<O.AF_unrel<<"\n";
	out<<"MDQ_rank="<<O.MDQ_rank<<"\n";
	out<<"MDQ="<<O.MDQ<<"\n";
	out<<"gene="<<O.gene<<"\n";
	
	
}


class calls_record{
	public:
	std::map<std::string,std::list<call_data> > L;
	
	void L_updater(const each_line_data& ELDO,const vcf_line_cols& VLCO,const parameters2& params,const double P_T,const double ExAC_AF_T,const double AF_unrel_T);
	double compare(const std::list<call_data>& C) const;
	
	double compare_aux(const std::vector<double>& C,const std::vector<double>& P) const;
};

double calls_record::compare_aux(const std::vector<double>& C,const std::vector<double>& P) const{
	// C is child, P is parent
	using namespace std;
	if( C.size()!=3 or P.size()!=3 )
		throw(domain_error("ERROR1 in calls_record::compare_aux"));
	
	vector<vector<double> > work_table(3,vector<double>(3,0));
	double max_val( -numeric_limits<double>::infinity() );
	for(int i(0);i<3;++i){
		for(int j(0);j<3;++j){
			work_table[i][j]+=( C[i]+P[j] );
			if( work_table[i][j]>max_val )
				max_val=work_table[i][j];
		}
	}
	double summ1(0);
	double summ2(0);
	for(int i(0);i<3;++i){
		for(int j(0);j<3;++j){
			work_table[i][j]-=max_val;
			work_table[i][j]=exp( work_table[i][j] );
			summ1+=work_table[i][j];
			if(j>=i)
				summ2+=work_table[i][j];
		}
	}
	if(summ2>summ1)
		throw(domain_error("ERROR2 in calls_record::compare_aux"));
	
	return (summ2/summ1);
}

double calls_record::compare(const std::list<call_data>& C) const{
	using namespace std;
	list<call_data>::const_iterator I(C.begin());
	double probty1(1);
	double probty2(1);
	for(;I!=C.end();++I){
		const call_data& CD( *I );
		const vector<double>& GL_L( CD.GL_L );
		const vector<double>& GL_L_P1( CD.GL_L_P1 );
		const vector<double>& GL_L_P2( CD.GL_L_P2 );
		probty1*= compare_aux(GL_L,GL_L_P1);
		probty2*= compare_aux(GL_L,GL_L_P2);
	}
	probty1=1-probty1;
	probty2=1-probty2;
	
	return (probty1*probty2);
}


void calls_record::L_updater(const each_line_data& ELDO,const vcf_line_cols& VLCO,const parameters2& param,const double P_T,const double ExAC_AF_T,const double AF_unrel_T){
	using namespace std;
	
	const double P( ELDO.HET_P_calc(VLCO) );
	//cout<<"P="<<P<<"\n";
	const double ExAC_AF( ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF );
	const double AF_unrel( param.AF_unrel );
	const int MDQ_rank( ELDO.INFO_data_obj.CSQ_data_obj.MDQ_rank );
	if( P>P_T && ExAC_AF<ExAC_AF_T && AF_unrel<AF_unrel_T && MDQ_rank>=17 ){
		const string gene( ELDO.INFO_data_obj.CSQ_data_obj.gene );
		call_data CD(VLCO,ELDO,param,P);
		L[gene].push_back(CD);
	}
}


std::ostream& operator<<(std::ostream& out,const calls_record& O){
	using namespace std;
	out.precision(30);
	map<string,list<call_data> >::const_iterator I( O.L.begin() );
	size_t count(0);
	for(;I!=O.L.end();++I){
		const list<call_data>& V( I->second );
		const double CC( O.compare(V) );
		if( V.size()>1 && CC>0.9 ){
			++count;
			out<<count<<") "<< I->first <<"\n";
			out<<"-----------\n";
			list<call_data>::const_iterator II( V.begin() );
			for(;II!=V.end();++II)
				out<< *II;
			out<<"\n";
		}
	}
	return out;
}


///////
//void AD_extractor(std::vector<std::vector<int> >& all_cols,const std::vector<std::string>& split_line,const size_t AD_col,const size_t ALT_count,const size_t FORMAT_length,const size_t FORMAT_col)




#endif
