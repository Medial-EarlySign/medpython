#if 0

// legacy code routines, kept - just in case....
//----------------------------------------------

//------------------------------------------------
void MedConvert::get_next_signal(vector<string> &buffered_lines, int &buffer_pos, ifstream &inf, int file_type, pid_data &curr, int &fpid, file_stat &curr_fstat, map<pair<string, string>, int>& missing_dict_vals)
{
	if (buffered_lines.empty()) {
		//fill buffer:
		if (read_file_to_buffer(inf, buffered_lines, read_lines_buffer)) {
			fpid = -1;
			--n_open_in_files;
			return;
		}
	}

	bool get_next = true;
	collected_data cd;
	GenericSigVec cd_sv;
	cd_sv.set_data(&cd.buf[0], 1);
	string vfield, vfield2;
	int i;
	int sid;
	char convert_mode = safe_mode > 0 ? 2 : 1;
	while (get_next) {
		if (buffer_pos >= buffered_lines.size()) {
			//need to refresh and read more/again from the buffer:
			buffer_pos = 0;
			buffered_lines.clear();
			if (read_lines_buffer > 0)
				buffered_lines.reserve(read_lines_buffer);
			bool finished = read_file_to_buffer(inf, buffered_lines, read_lines_buffer);
			if (finished) {
				get_next = false;
				fpid = -1;
				--n_open_in_files;
				return;
			}
		}

		if (buffer_pos < buffered_lines.size()) {
			string &curr_line = buffered_lines[buffer_pos];
			++buffer_pos;
			//remove first row
			/*if (curr_line[curr_line.size() - 1] == '\r')
			curr_line.erase(curr_line.size() - 1);
			boost::trim(curr_line);*/
			curr_fstat.n_lines++;
			//if ((curr_line.size() > 1) && (curr_line[0] != '#')) {
			//if (fpid == 5025392)
			//	MLOG("fpid %d : %s\n", fpid, curr_line.c_str());
			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));
			if (fields.size() == 1)
				split(fields, curr_line, boost::is_any_of(" ")); // bypass for files that were space delimited

			curr_fstat.n_relevant_lines++;
			//if (fpid == 5025392)
			//				MLOG("working on: (fpid %d) (curr.pid %d) (file_type %d) (f[0] %s) (nfields %d) ##>%s<##\n",fpid,curr.pid,file_type,fields[0].c_str(),fields.size(),curr_line.c_str());
			//if (fields.size() > 6) MLOG("WEIRD f[6]= ##>%s<##\n",fields[5].c_str());
			if (((file_type == 1) && (fields.size() == 4)) ||
				((file_type == 2) && (fields.size() >= 3)) ||
				((file_type == 3) && (fields.size() >= 3))) {

				int line_pid;
				try {
					line_pid = med_stoi(fields[0]);
				}
				catch (...) {
					MERR("ERROR: bad format in file %s with first token of pid, in line %d:\n%s\n",
						curr_fstat.fname.c_str(), curr_fstat.n_parsed_lines, curr_line.c_str());
					throw;
				}
				//if (fpid == 5025392)
				//	MLOG("working on: (fpid %d) (curr.pid %d) (file_type %d) (line_pid %d) %s\n",fpid,curr.pid,file_type,line_pid,curr_line.c_str());
				if (line_pid == curr.pid) {
					cd.zero();
					if (file_type == 1) {

						// Registry file //format: pid , stage(string) , date, location(number)	// tab delimited 

						try {
							cd_sv.init_from_sigtype(T_DateVal);
							// Cancer_Location
							i = sid2serial[dict.id(string("Cancer_Location"))];
							cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(default_time_unit, fields[2], convert_mode));
							cd_sv.setVal(0, 0, (float)(dict.id(fields[3])));
							//#pragma omp critical
							curr.raw_data[i].push_back(cd);

							// Cancer_Stage
							i = sid2serial[dict.id(string("Cancer_Stage"))];
							cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(default_time_unit, fields[2], convert_mode));
							cd_sv.setVal(0, 0, (float)(med_stoi(fields[1])));
							//#pragma omp critical
							curr.raw_data[i].push_back(cd);

							++curr_fstat.n_parsed_lines;
						}
						catch (...) {
							MERR("ERROR: bad format in parsing registry file %s in line %d:\n%s\n",
								curr_fstat.fname.c_str(), curr_fstat.n_parsed_lines, curr_line.c_str());
							throw;
						}

					}
					else if (file_type >= 2) {
						if (codes2names.find(fields[1]) == codes2names.end())
							MTHROW_AND_ERR("MedConvert: ERROR: unrecognized signal name %s (need to add to codes_to_signals file) in file %s :: curr_line is %s\n",
								fields[1].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
						sid = dict.id(codes2names[fields[1]]);
						if (sid < 0)
							MTHROW_AND_ERR("MedConvert: ERROR: signal name %s converted to %s is not in dict in file %s :: curr_line is %s\n",
								fields[1].c_str(), codes2names[fields[1]].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
						//MLOG("here001 %s %d %d \n", codes2names[fields[1]].c_str(), sid, sids_to_load[sid]);
						if (!sids_to_load[sid])
							continue;
						int section = dict.section_id(sigs.name(sid));
						try {
							i = sid2serial[sid];
							SignalInfo& info = sigs.Sid2Info[sid];
							cd_sv.init(info);
							if (cd_sv.size() > MAX_COLLECTED_DATA_SIZE) {
								MTHROW_AND_ERR("ERROR: cd_sv.size() (%d) > MAX_COLLECTED_DATA_SIZE (%d), Please Increase MAX_COLLECTED_DATA_SIZE\n", (int)cd_sv.size(), (int)MAX_COLLECTED_DATA_SIZE);
							}
							int time_unit = info.time_unit == MedTime::Undefined ? default_time_unit : info.time_unit;
							if (file_type == 3) {
								// backward compatibility - if file type is DATA_S (3), then all val channels are assumed categorical
								for (int j = 0; j < info.n_val_channels; j++) {
									if (info.is_categorical_per_val_channel[j] != 1) {
#pragma omp critical
										info.is_categorical_per_val_channel[j] = 1;
									}
								}
							}
							switch (sigs.type(sid)) {

							case T_Value:
								//cd.date = 0;
								if (fields.size() == 3) {
									if (sigs.is_categorical_channel(sid, 0))
										cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[2]));
									else cd_sv.setVal(0, 0, med_stof(fields[2]));
								}
								else { // backward compatible with date 0 trick to load value only data
									if (sigs.is_categorical_channel(sid, 0))
										cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[3]));
									else
										cd_sv.setVal(0, 0, med_stof(fields[3]));
								}
								break;

							case T_DateVal:
								cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal(0, 0, med_stof(fields[3]));
								break;

							case T_TimeVal:
								cd_sv.setTime<long long>(0, 0, stoll(fields[2]));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal(0, 0, med_stof(fields[3]));
								break;

							case T_DateRangeVal:
								cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								cd_sv.setTime(0, 1, med_time_converter.convert_datetime_safe(time_unit, fields[3], convert_mode));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal(0, 0, med_stof(fields[4]));
								break;

							case T_TimeStamp:
								cd_sv.setTime<long long>(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								break;

							case T_TimeRangeVal:
								cd_sv.setTime<long long>(0, 0, stoll(fields[2]));
								cd_sv.setTime<long long>(0, 1, stoll(fields[3]));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal(0, 0, med_stof(fields[4]));
								break;

							case T_DateVal2:
								cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal(0, 0, med_stof(fields[3]));
								if (sigs.is_categorical_channel(sid, 1))
									cd_sv.setVal<unsigned short>(0, 1, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal<unsigned short>(0, 1, (unsigned short)med_stoi(fields[4]));
								break;

							case T_TimeLongVal:

								cd_sv.setTime<long long>(0, 0, stoll(fields[2]));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal<long long>(0, 0, (long long)dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal<long long>(0, 0, med_stof(fields[3]));
								break;

							case T_DateShort2:
								cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal<short>(0, 0, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal<short>(0, 0, med_stof(fields[3]));
								if (sigs.is_categorical_channel(sid, 1))
									cd_sv.setVal<short>(0, 1, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal<short>(0, 1, med_stof(fields[4]));
								break;

							case T_ValShort2:
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal<short>(0, 0, dict.get_id_or_throw(section, fields[2]));
								else cd_sv.setVal<short>(0, 0, med_stof(fields[2]));
								if (sigs.is_categorical_channel(sid, 1))
									cd_sv.setVal<short>(0, 1, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal<short>(0, 1, med_stof(fields[3]));
								break;

							case T_ValShort4:
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal<short>(0, 0, dict.get_id_or_throw(section, fields[2]));
								else cd_sv.setVal<short>(0, 0, med_stof(fields[2]));
								if (sigs.is_categorical_channel(sid, 1))
									cd_sv.setVal<short>(0, 1, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal<short>(0, 1, med_stof(fields[3]));
								if (sigs.is_categorical_channel(sid, 2))
									cd_sv.setVal<short>(0, 2, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal<short>(0, 2, med_stof(fields[4]));
								if (sigs.is_categorical_channel(sid, 3))
									cd_sv.setVal<short>(0, 3, dict.get_id_or_throw(section, fields[5]));
								else cd_sv.setVal<short>(0, 3, med_stof(fields[5]));
								break;
							case T_CompactDateVal:
								cd_sv.setTime<unsigned short>(0, 0, (int)med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal<unsigned short>(0, 0, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal<unsigned short>(0, 0, (unsigned short)med_stoi(fields[3]));
								break;

							case T_DateRangeVal2:
								cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								cd_sv.setTime(0, 1, med_time_converter.convert_datetime_safe(time_unit, fields[3], convert_mode));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal(0, 0, med_stof(fields[4]));
								if (sigs.is_categorical_channel(sid, 1))
									cd_sv.setVal(0, 1, dict.get_id_or_throw(section, fields[5]));
								else cd_sv.setVal(0, 1, med_stof(fields[5]));
								break;

							case T_DateFloat2:
								cd_sv.setTime(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));
								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal(0, 0, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal(0, 0, med_stof(fields[3]));
								if (sigs.is_categorical_channel(sid, 1))
									cd_sv.setVal(0, 1, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal(0, 1, med_stof(fields[4]));
								break;

							case T_TimeShort4:

								cd_sv.setTime<long long>(0, 0, med_time_converter.convert_datetime_safe(time_unit, fields[2], convert_mode));

								if (sigs.is_categorical_channel(sid, 0))
									cd_sv.setVal<unsigned short>(0, 0, dict.get_id_or_throw(section, fields[3]));
								else cd_sv.setVal<unsigned short>(0, 0, med_stof(fields[3]));

								if (sigs.is_categorical_channel(sid, 1))
									cd_sv.setVal<unsigned short>(0, 1, dict.get_id_or_throw(section, fields[4]));
								else cd_sv.setVal<unsigned short>(0, 1, med_stof(fields[4]));

								if (sigs.is_categorical_channel(sid, 2))
									cd_sv.setVal<unsigned short>(0, 2, dict.get_id_or_throw(section, fields[5]));
								else cd_sv.setVal<unsigned short>(0, 2, med_stof(fields[5]));

								if (sigs.is_categorical_channel(sid, 3))
									cd_sv.setVal<unsigned short>(0, 3, dict.get_id_or_throw(section, fields[6]));
								else cd_sv.setVal<unsigned short>(0, 3, med_stof(fields[6]));

								break;
							case T_Generic:
							{
								int field_i = 2;
								for (int tchan = 0; tchan < cd_sv.n_time; tchan++) {
									switch (cd_sv.time_channel_types[tchan]) {

									case GenericSigVec::type_enc::UINT8:   //unsigned char
									case GenericSigVec::type_enc::UINT32:  //unsigned int
									case GenericSigVec::type_enc::UINT64:  //unsigned long long
									case GenericSigVec::type_enc::INT8:    //char
									case GenericSigVec::type_enc::INT16:   //short
									case GenericSigVec::type_enc::FLOAT32: //float
									case GenericSigVec::type_enc::FLOAT64: //double
									case GenericSigVec::type_enc::FLOAT80: //long double

									case GenericSigVec::type_enc::INT32:   //int
										cd_sv.setTime(0, tchan, med_time_converter.convert_datetime_safe(time_unit, fields[field_i], convert_mode));
										break;
										//TODO: figure out when to use stoll and time_converter
									case GenericSigVec::type_enc::INT64:   //long long
																		   //cd_sv.setTime<long long>(0, tchan, stoll(fields[field_i]));
										cd_sv.setTime<long long>(0, tchan, med_time_converter.convert_datetime_safe(time_unit, fields[field_i], convert_mode));
										break;
									case GenericSigVec::type_enc::UINT16:  //unsigned short
										if (1) {
											int value = (int)med_time_converter.convert_datetime_safe(time_unit, fields[field_i], convert_mode);
											if (value < 0)
												MTHROW_AND_ERR("MedConvert: get_next_signal: Detected attempt to assign negative number (%d) into unsigned time channel %d :: curr_line is '%s'\n", value, tchan, curr_line.c_str());
											cd_sv.setTime<unsigned short>(0, tchan, value);
										}
										break;
									}
									field_i++;
								}
								for (int vchan = 0; vchan < cd_sv.n_val; vchan++) {
									switch (cd_sv.val_channel_types[vchan]) {

									case GenericSigVec::type_enc::UINT8:   //unsigned char
									case GenericSigVec::type_enc::UINT32:  //unsigned int
									case GenericSigVec::type_enc::INT8:    //char
									case GenericSigVec::type_enc::INT32:   //int
									case GenericSigVec::type_enc::INT64:   //long long
									case GenericSigVec::type_enc::FLOAT64: //double
									case GenericSigVec::type_enc::FLOAT80: //long double

									case GenericSigVec::type_enc::FLOAT32: //float
										if (sigs.is_categorical_channel(sid, vchan))
											cd_sv.setVal(0, vchan, dict.get_id_or_throw(section, fields[field_i]));
										else cd_sv.setVal(0, vchan, med_stof(fields[field_i]));
										break;

									case GenericSigVec::type_enc::UINT16:  //unsigned short
										if (sigs.is_categorical_channel(sid, vchan))
											cd_sv.setVal<unsigned short>(0, vchan, dict.get_id_or_throw(section, fields[field_i]));
										else {
											auto value = med_stoi(fields[field_i]);
											if (value < 0)
												MTHROW_AND_ERR("MedConvert: get_next_signal: Detected attempt to assign negative number (%d) into unsigned value channel %d :: curr_line is '%s'\n", value, vchan, curr_line.c_str());
											cd_sv.setVal<unsigned short>(0, vchan, (unsigned short)value);
										}
										break;
									case GenericSigVec::type_enc::UINT64:  //unsigned long long
										if (sigs.is_categorical_channel(sid, vchan))
											cd_sv.setVal<long long>(0, vchan, (long long)dict.get_id_or_throw(section, fields[field_i]));
										else cd_sv.setVal<long long>(0, vchan, med_stof(fields[field_i]));
										break;
									case GenericSigVec::type_enc::INT16:   //short
										if (sigs.is_categorical_channel(sid, vchan))
											cd_sv.setVal<short>(0, vchan, dict.get_id_or_throw(section, fields[field_i]));
										else cd_sv.setVal<short>(0, vchan, med_stof(fields[field_i]));
										break;
									}
									field_i++;
								}

								break;
							}
							default:
								MTHROW_AND_ERR("MedConvert: get_next_signal: unknown signal type %d for sid %d\n",
									sigs.type(sid), sid);
							}


							//#pragma omp critical
							curr.raw_data[i].push_back(cd);
							curr_fstat.n_parsed_lines++;
						}
						catch (invalid_argument e) {
							pair<string, string> my_key = make_pair(sigs.name(sid), string(e.what()));
							if (missing_dict_vals.find(my_key) == missing_dict_vals.end()) {
#pragma omp critical
								missing_dict_vals[my_key] = 1;
							}
							else {
#pragma omp atomic
								++missing_dict_vals[my_key];
							}
							if (missing_dict_vals.size() < 10)
								MWARN("MedConvert::get_next_signal: missing from dictionary (sig [%s], type %d) : file [%s] : line [%s] \n",
									sigs.name(sid).c_str(), sigs.type(sid), curr_fstat.fname.c_str(), curr_line.c_str());
							if (!full_error_file.empty()) {
								if (err_log_file.good()) {
									err_log_file << "MISSING_FROM_DICTIONARY"
										<< "\t" << sigs.name(sid) << " (" << sigs.type(sid) << ")"
										<< "\t" << curr_fstat.fname
										<< "\t" << "\"" << curr_line << "\""
										<< "\n";
								}
							}
						}
						catch (...) {
							curr_fstat.n_bad_format_lines++;
							if (curr_fstat.n_bad_format_lines < 10) {
								MWARN("MedConvert::get_next_signal: bad format in parsing file %s (file_type=%d) in line %d:\n%s\n",
									curr_fstat.fname.c_str(), file_type, curr_fstat.n_parsed_lines, curr_line.c_str());
							}
							if (!full_error_file.empty()) {
								if (err_log_file.good()) {
									err_log_file << "BAD_FILE_FORMAT"
										<< "\t" << curr_fstat.fname << " (" << file_type << ")"
										<< "\t" << curr_fstat.n_parsed_lines
										<< "\t" << "\"" << curr_line << "\""
										<< "\n";
								}
							}
						}
					}
				}
				else if (line_pid < fpid) {
					MWARN("MedConvert: get_next_signal: fpid is %d , but got line: %s\n", fpid, curr_line.c_str());
					if (safe_mode)
						MTHROW_AND_ERR("MedConvert: ERROR: file %s seems to be not sorted by pid\n", curr_fstat.fname.c_str());
				}
				else {
					fpid = line_pid;
					--buffer_pos; // roll file back to the start of curr line
								  //vector<string> sub_lines(buffered_lines.begin() + buffer_pos, buffered_lines.end());
								  //buffered_lines.swap(sub_lines);
								  //update file stats
					--curr_fstat.n_lines;
					--curr_fstat.n_relevant_lines;
					get_next = false;
				}
			}
			//}
		}
		else
			get_next = false;
	}
	if (fpid > MAX_PID_TO_TAKE) {
		fpid = -1;
		--n_open_in_files;

	}
}


//------------------------------------------------
void MedConvert::get_next_signal_new_modes(vector<string> &buffered_lines, int &buffer_pos, ifstream &inf, int file_type, pid_data &curr, int &fpid, file_stat &curr_fstat, map<pair<string, string>, int>& missing_dict_vals)
{
	bool get_next = true;
	collected_data cd;
	GenericSigVec cd_sv;
	cd_sv.set_data(&cd.buf[0], 1);
	string vfield, vfield2;
	int i;
	int sid;

	while (get_next) {
		if (buffer_pos >= buffered_lines.size() || buffered_lines.empty()) {
			//need to refresh and read more/again from the buffer:
			buffer_pos = 0;
			buffered_lines.clear();
			if (read_file_to_buffer(inf, buffered_lines, read_lines_buffer)) {
				get_next = false;
				fpid = -1;
				--n_open_in_files;
				return;
			}
		}

		if (buffer_pos < buffered_lines.size()) {
			string &curr_line = buffered_lines[buffer_pos];
			++buffer_pos;
			curr_fstat.n_lines++;
			vector<string> fields;
			split(fields, curr_line, boost::is_any_of("\t"));

			curr_fstat.n_relevant_lines++;
			int line_pid;
			try {
				line_pid = med_stoi(fields[0]);
			}
			catch (...) {
				MERR("ERROR: bad format in file %s with first token of pid, in line %d:\n%s\n",
					curr_fstat.fname.c_str(), curr_fstat.n_parsed_lines, curr_line.c_str());
				throw;
			}

			if (line_pid == curr.pid) {
				cd.zero();
				if (codes2names.find(fields[1]) == codes2names.end())
					MTHROW_AND_ERR("MedConvert: ERROR: unrecognized signal name %s (need to add to codes_to_signals file) in file %s :: curr_line is %s\n",
						fields[1].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());
				sid = dict.id(codes2names[fields[1]]);
				if (sid < 0)
					MTHROW_AND_ERR("MedConvert: ERROR: signal name %s converted to %s is not in dict in file %s :: curr_line is %s\n",
						fields[1].c_str(), codes2names[fields[1]].c_str(), curr_fstat.fname.c_str(), curr_line.c_str());

				if (!sids_to_load[sid])
					continue;

				try {
					i = sid2serial[sid];
					parse_fields_into_gsv(curr_line, fields, sid, cd_sv);
#pragma omp critical
					{
						curr.raw_data[i].push_back(cd);
						curr_fstat.n_parsed_lines++;
					}
				}
				catch (invalid_argument e) {
					pair<string, string> my_key = make_pair(sigs.name(sid), string(e.what()));
					if (missing_dict_vals.find(my_key) == missing_dict_vals.end()) {
#pragma omp critical
						missing_dict_vals[my_key] = 1;
					}
					else {
#pragma omp atomic
						++missing_dict_vals[my_key];
					}
					if (missing_dict_vals.size() < 10)
						MWARN("MedConvert::get_next_signal: missing from dictionary (sig [%s], type %d) : file [%s] : line [%s] \n",
							sigs.name(sid).c_str(), sigs.type(sid), curr_fstat.fname.c_str(), curr_line.c_str());
					if (!full_error_file.empty()) {
						if (err_log_file.good()) {
							err_log_file << "MISSING_FROM_DICTIONARY"
								<< "\t" << sigs.name(sid) << " (" << sigs.type(sid) << ")"
								<< "\t" << curr_fstat.fname
								<< "\t" << "\"" << curr_line << "\""
								<< "\n";
						}
					}
				}
				catch (...) {
					curr_fstat.n_bad_format_lines++;
					if (curr_fstat.n_bad_format_lines < 10) {
						MWARN("MedConvert::get_next_signal: bad format in parsing file %s (file_type=%d) in line %d:\n%s\n",
							curr_fstat.fname.c_str(), file_type, curr_fstat.n_parsed_lines, curr_line.c_str());
					}
					if (!full_error_file.empty()) {
						if (err_log_file.good()) {
							err_log_file << "BAD_FILE_FORMAT"
								<< "\t" << curr_fstat.fname << " (" << file_type << ")"
								<< "\t" << curr_fstat.n_parsed_lines
								<< "\t" << "\"" << curr_line << "\""
								<< "\n";
						}
					}
				}
			}
			else if (line_pid < fpid) {
				MWARN("MedConvert: get_next_signal: fpid is %d , but got line: %s\n", fpid, curr_line.c_str());
				if (safe_mode)
					MTHROW_AND_ERR("MedConvert: ERROR: file %s seems to be not sorted by pid\n", curr_fstat.fname.c_str());
			}
			else {
				fpid = line_pid;
				--buffer_pos; // roll file back to the start of curr line
				--curr_fstat.n_lines;
				--curr_fstat.n_relevant_lines;
				get_next = false;
			}
		}
		else
			get_next = false;
	}

	if (fpid > MAX_PID_TO_TAKE) {
		fpid = -1;
		--n_open_in_files;

	}
}

//------------------------------------------------
void merge_changes(const pid_data &curr_i, pid_data &curr) {
#pragma omp critical
	for (size_t i = 0; i < curr.raw_data.size(); ++i) {
		if (!curr_i.raw_data[i].empty())
			curr.raw_data[i].insert(curr.raw_data[i].end(), curr_i.raw_data[i].begin(), curr_i.raw_data[i].end());
	}
}

//------------------------------------------------
int MedConvert::write_indexes(pid_data &curr)
{
	if (curr.pid < 0)
		MTHROW_AND_ERR("MedConvert::write_indexes negative pid %d", curr.pid);
	// first we sort all elements by time
	for (int i = 0; i < curr.raw_data.size(); i++) {
		GenericSigVec gsv1;
		auto& info = sigs.Sid2Info[serial2siginfo[i].sid];
		gsv1.init(info);
		sort(curr.raw_data[i].begin(), curr.raw_data[i].end(),
			[&](const collected_data &v1, const collected_data &v2)
		{
			return gsv1.compareTimeLt(&v1.buf[0], 0, &v2.buf[0], 0);
		});
	}

	// getting rid of duplicates
	vector<collected_data>::iterator it;
	for (int i = 0; i < curr.raw_data.size(); i++) {
		int struct_size = sigs.Sid2Info[serial2siginfo[i].sid].bytes_len;
		it = unique(curr.raw_data[i].begin(), curr.raw_data[i].end(), [=](const collected_data &v1, const collected_data &v2) {
			return memcmp(v1.buf, v2.buf, struct_size) == 0;
			//return cd_to_tuple(v1) == cd_to_tuple(v2); 
		});
		curr.raw_data[i].resize(distance(curr.raw_data[i].begin(), it));
	}


	// sanity checks - things we force to have, and things we force to have as single.
	//	if (curr.raw_data[sid2serial[dict.id(string("GENDER"))]].size() != 1) return -1;
	//	if (curr.raw_data[sid2serial[dict.id(string("BYEAR"))]].size() != 1) return -1;

	// forced signals
	for (int i = 0; i < forced.size(); i++) {
		if (curr.raw_data[sid2serial[dict.id(forced[i])]].size() != 1) {
			if (missing_forced_signals.find(forced[i]) == missing_forced_signals.end())
				missing_forced_signals[forced[i]] = 1;
			else
				missing_forced_signals[forced[i]] += 1;
			if (missing_forced_signals[forced[i]] < 10)
				MLOG("MedConvert: pid %d is missing forced signal %s (%d,%d,%d)\n", curr.pid, forced[i].c_str(),
					dict.id(forced[i]), sid2serial[dict.id(forced[i])], curr.raw_data[sid2serial[dict.id(forced[i])]].size());

			if (!full_error_file.empty()) {
				if (err_log_file.good()) {
					err_log_file << "MISSING_FORCED_SIGNAL"
						<< "\t" << curr.pid
						<< "\t" << forced[i] << "(" << dict.id(forced[i]) << ","
						<< sid2serial[dict.id(forced[i])] << ","
						<< curr.raw_data[sid2serial[dict.id(forced[i])]].size() << ")"
						<< "\n";
				}
			}
			return -1;
		}
	}

	// writing indexes
	int fno;
	int n_pid_sigs;
	if (test_run_max_pids == 0)
		for (fno = 0; fno < index_fnames.size(); fno++)
			if (data_f[fno] != NULL)
			{
				n_pid_sigs = 0;
				if (mode < 3) {
					for (int i = 0; i < curr.raw_data.size(); i++)
						if (curr.raw_data[i].size() > 0 && serial2siginfo[i].fno == fno &&
							(serial2siginfo[i].type >= 0 && serial2siginfo[i].type < T_Last))
							//if (curr.raw_data[i].size() > 0 && sid2fno[serial2sid[i]] == fno &&
							//(sigs.type(serial2sid[i]) >= 0 &&   sigs.type(serial2sid[i])<T_Last))
							n_pid_sigs++;
				}
				else {
					// in this mode fno is i... and there's one option for it
					if (curr.raw_data[fno].size() > 0 && serial2siginfo[fno].type >= 0 && serial2siginfo[fno].type < T_Last)
						n_pid_sigs++;
					//MLOG("i=%d/%d n_pid_sigs %d\n", fno, index_fnames.size(), n_pid_sigs);
				}


				if (n_pid_sigs > 0) {

					// write packet header: magic number + pid + number of signals
					unsigned long long magic = MED_MAGIC_NUM;
					int pid = curr.pid;
					if (mode < 3) {
						index_f[fno]->write((char *)&magic, sizeof(unsigned long long));
						index_f[fno]->write((char *)&pid, sizeof(int));
						index_f[fno]->write((char *)&n_pid_sigs, sizeof(int));
					}

					// write data and index pointer for each signal
					for (int i = 0; i < curr.raw_data.size(); i++) {

						int ilen = (int)curr.raw_data[i].size();
						if (ilen > 0) {
							int sid = serial2sid[i];
							int sid_type = serial2siginfo[i].type; //sigs.type(sid);
							int sid_fno = serial2siginfo[i].fno; //sid2fno[sid];

							if ((ilen > 0) && (sid_fno == fno) && (sid_type >= 0 && sid_type < T_Last)) {

								//int sid = serial2sid[i];
								unsigned short file_n = fno;
								unsigned long long pos = data_f_pos[fno];
								int len = 0;

								int struct_len = (int)sigs.Sid2Info[serial2siginfo[i].sid].bytes_len;
								len = struct_len * ilen;
								for (int j = 0; j < ilen; j++) {
									data_f[fno]->write((char *)&(curr.raw_data[i][j].buf[0]), struct_len);
								}


								//MLOG("writing to fno %d : sid %d file_n %d pos %ld len %d\n", fno, sid, file_n, pos, len);
								if (mode < 3) {
									index_f[fno]->write((char *)&sid, sizeof(int));
									index_f[fno]->write((char *)&file_n, sizeof(short));
									index_f[fno]->write((char *)&pos, sizeof(unsigned long long));
									index_f[fno]->write((char *)&len, sizeof(int));
								}
								else {
									indexes[fno].insert(pid, ilen);
								}
								data_f_pos[fno] += len;

							}
						}
					}
				}

			}
	return 0;
}



// code from create_indexes:....
//-------------------------------
/*


#pragma omp parallel for schedule(dynamic,dynamic_threads)
for (int i = 0; i < n_files_opened; i++) {
int fpid = c_pid;
curr_f[i].pid = c_pid;
curr_f[i].raw_data.clear();
curr_f[i].raw_data.resize(serial2sid.size());
if (pid_in_file[i] <= c_pid) {
get_next_signal_new_modes(file_to_lines[i], file_buffer_pos[i], infs[i], file_type[i], curr_f[i], fpid, fstats[i], missing_dict_vals);
pid_in_file[i] = fpid; // current pid after the one we wanted
}
}


#pragma omp parallel for schedule(dynamic,dynamic_threads)
for (int j = 0; j < curr.raw_data.size(); j++)
for (int i = 0; i < n_files_opened; i++)
if (curr_f[i].raw_data[j].size() > 0)
curr.raw_data[j].insert(curr.raw_data[j].end(), curr_f[i].raw_data[j].begin(), curr_f[i].raw_data[j].end());
*/

/*
//bool mark_err = false;
#pragma omp parallel for schedule(dynamic) if (run_parallel)
if (pid_in_file[i] <= c_pid) {
if (run_parallel) {
pid_data curr_i;
curr_i.pid = c_pid;
curr_i.raw_data.resize(serial2sid.size());
//MLOG("file %d :: pid_int_file %d fpid %d\n", i, pid_in_file[i], fpid);
get_next_signal(file_to_lines[i], file_buffer_pos[i], infs[i], file_type[i], curr_i, fpid, fstats[i], missing_dict_vals);
pid_in_file[i] = fpid; // current pid after the one we wanted
//merge into curr from curr_i:collect_lines
merge_changes(curr_i, curr);
}
else {
//update pid if exit code is zero
get_next_signal(file_to_lines[i], file_buffer_pos[i], infs[i], file_type[i], curr, fpid, fstats[i], missing_dict_vals);
pid_in_file[i] = fpid; // current pid after the one we wanted
}

}
}
//MLOG("i=%d c_pid=%d fpid=%d curr %d %d %d\n",i,c_pid,fpid,curr.pid,n_files_opened,n_open_in_files);
*/

// code from write_indexes:....
//------------------------------
/*
if (sid_type == T_Value) {
len = (int)sizeof(SVal)*ilen;
SVal sv;
for (int j = 0; j < ilen; j++) {
sv.val = curr.raw_data[i][j].val;
data_f[fno]->write((char *)&sv, sizeof(SVal));
}
}

if (sid_type == T_DateVal) {
len = (int)sizeof(SDateVal)*ilen;
SDateVal sdv;
for (int j = 0; j < ilen; j++) {
sdv.date = curr.raw_data[i][j].date;
sdv.val = curr.raw_data[i][j].val;
data_f[fno]->write((char *)&sdv, sizeof(SDateVal));
}
}

if (sid_type == T_DateRangeVal) {
len = (int)sizeof(SDateRangeVal)*ilen;
SDateRangeVal sdrv;
for (int j = 0; j < ilen; j++) {
sdrv.date_start = curr.raw_data[i][j].date;
sdrv.date_end = curr.raw_data[i][j].date2;
sdrv.val = curr.raw_data[i][j].val;
data_f[fno]->write((char *)&sdrv, sizeof(SDateRangeVal));
}
}

if (sid_type == T_DateRangeVal2) {
len = (int)sizeof(SDateRangeVal2)*ilen;
SDateRangeVal2 sdrv;
for (int j = 0; j < ilen; j++) {
sdrv.date_start = curr.raw_data[i][j].date;
sdrv.date_end = curr.raw_data[i][j].date2;
sdrv.val = curr.raw_data[i][j].val;
sdrv.val2 = curr.raw_data[i][j].f_val2;
data_f[fno]->write((char *)&sdrv, sizeof(SDateRangeVal2));
}
}

if (sid_type == T_DateFloat2) {
len = (int)sizeof(SDateFloat2)*ilen;
SDateFloat2 sdrv;
for (int j = 0; j < ilen; j++) {
sdrv.date = curr.raw_data[i][j].date;
sdrv.val = curr.raw_data[i][j].val;
sdrv.val2 = curr.raw_data[i][j].f_val2;
data_f[fno]->write((char *)&sdrv, sizeof(SDateFloat2));
}
}

if (sid_type == T_TimeVal) {
len = (int)sizeof(STimeVal)*ilen;
STimeVal stv;
for (int j = 0; j < ilen; j++) {
stv.time = curr.raw_data[i][j].time;
stv.val = curr.raw_data[i][j].val;
data_f[fno]->write((char *)&stv, sizeof(STimeVal));
}
}

if (sid_type == T_TimeRangeVal) {
len = (int)sizeof(STimeRangeVal)*ilen;
STimeRangeVal strv;
for (int j = 0; j < ilen; j++) {
strv.time_start = curr.raw_data[i][j].time;
strv.time_end = curr.raw_data[i][j].time2;
strv.val = curr.raw_data[i][j].val;
data_f[fno]->write((char *)&strv, sizeof(STimeRangeVal));
}
}

if (sid_type == T_TimeStamp) {
len = (int)sizeof(STimeStamp)*ilen;
STimeStamp sts;
for (int j = 0; j < ilen; j++) {
sts.time = curr.raw_data[i][j].time;
data_f[fno]->write((char *)&sts, sizeof(STimeStamp));
}
}

if (sid_type == T_DateVal2) {
len = (int)sizeof(SDateVal2)*ilen;
SDateVal2 sdv2;
for (int j = 0; j < ilen; j++) {
sdv2.date = curr.raw_data[i][j].date;
sdv2.val = curr.raw_data[i][j].val;
sdv2.val2 = curr.raw_data[i][j].val2;
data_f[fno]->write((char *)&sdv2, sizeof(SDateVal2));
}
}

if (sid_type == T_TimeLongVal) {
len = (int)sizeof(STimeLongVal)*ilen;
STimeLongVal stv;
for (int j = 0; j < ilen; j++) {
stv.time = curr.raw_data[i][j].time;
stv.val = curr.raw_data[i][j].longVal;
data_f[fno]->write((char *)&stv, sizeof(STimeLongVal));
}
}

if (sid_type == T_DateShort2) {
len = (int)sizeof(SDateShort2)*ilen;
SDateShort2 sds2;
for (int j = 0; j < ilen; j++) {
sds2.date = curr.raw_data[i][j].date;
sds2.val1 = curr.raw_data[i][j].val1;
sds2.val2 = curr.raw_data[i][j].val2;
data_f[fno]->write((char *)&sds2, sizeof(SDateShort2));
}
}

if (sid_type == T_ValShort2) {
len = (int)sizeof(SValShort2)*ilen;
SValShort2 svs2;
for (int j = 0; j < ilen; j++) {
svs2.val1 = curr.raw_data[i][j].val1;
svs2.val2 = curr.raw_data[i][j].val2;
data_f[fno]->write((char *)&svs2, sizeof(SValShort2));
}
}

if (sid_type == T_ValShort4) {
len = (int)sizeof(SValShort4)*ilen;
SValShort4 svs4;
for (int j = 0; j < ilen; j++) {
svs4.val1 = curr.raw_data[i][j].val1;
svs4.val2 = curr.raw_data[i][j].val2;
svs4.val3 = curr.raw_data[i][j].val3;
svs4.val4 = curr.raw_data[i][j].val4;
data_f[fno]->write((char *)&svs4, sizeof(SValShort4));
}
}

if (sid_type == T_CompactDateVal) {
len = (int)sizeof(SCompactDateVal)*ilen;
SCompactDateVal scdv;
for (int j = 0; j < ilen; j++) {
scdv.compact_date = date_to_compact_date(curr.raw_data[i][j].date);
scdv.val = curr.raw_data[i][j].val1;
data_f[fno]->write((char *)&scdv, sizeof(SCompactDateVal));
}
}

if (sid_type == T_TimeShort4) {
len = (int)sizeof(STimeShort4)*ilen;
STimeShort4 sts4;
for (int j = 0; j < ilen; j++) {
sts4.time = curr.raw_data[i][j].time;
sts4.val1 = curr.raw_data[i][j].val1;
sts4.val2 = curr.raw_data[i][j].val2;
sts4.val3 = curr.raw_data[i][j].val3;
sts4.val4 = curr.raw_data[i][j].val4;
data_f[fno]->write((char *)&sts4, sizeof(STimeShort4));
}
}
*/



#endif