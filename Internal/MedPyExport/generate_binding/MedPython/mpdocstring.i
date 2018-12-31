%define MPDOCSTRING
"
    class Dictionary
        Methods:
            __init__(PidRepository rep) -> Dictionary
            get_members_to_all_sets(int section_id, int * members_array) -> IntVecIntMapAdaptor
            id(string & signame) -> int
            name(int id) -> string
            prep_sets_lookup_table(int section_id, StringVector set_names)
            section_id(string name) -> int


    class FeatureAttr
        Methods:
            imputed -> bool   (property getter)
            normalized -> bool   (property getter)
            imputed <- bool   (property setter)
            normalized <- bool   (property setter)
            __init__() -> FeatureAttr

        Data descriptors:
            imputed -> bool   (property getter)
            normalized -> bool   (property getter)


    class Features
        Methods:
            attributes -> StringFeatureAttrMapAdaptor   (property getter)
            data -> StringVecFloatMapAdaptor   (property getter)
            pid_pos_len -> IntPairIntIntMapAdaptor   (property getter)
            samples -> SampleVectorAdaptor   (property getter)
            tags -> StringUOSetStringMapAdaptor   (property getter)
            time_unit -> int   (property getter)
            weights    (property getter)
            time_unit <- int   (property setter)
            __init__() -> Features
            __init__(int _time_unit) -> Features
            __init__(Features other) -> Features
            append_samples(IdSamples in_samples)
            append_samples(Samples in_samples)
            clear()
            filter(StringVector selectedFeatures) -> int
            get_as_matrix(Mat mat)
            get_as_matrix(Mat mat, vector< string > names)
            get_as_matrix(Mat mat, vector< string > const names, int * int_in_buf)
            get_crc() -> unsigned int
            get_feature_names() -> StringVector
            get_max_serial_id_cnt() -> int
            get_pid_len(int pid) -> int
            get_pid_pos(int pid) -> int
            get_samples(Samples outSamples)
            init_all_samples(IdSamplesVectorAdaptor in_samples)
            init_pid_pos_len()
            insert_samples(IdSamples in_samples, int index)
            print_csv()
            read_from_csv_mat(string const & csv_fname) -> int
            set_as_matrix(Mat mat)
            set_time_unit(int _time_unit)
            version() -> int
            write_as_csv_mat(string const & csv_fname) -> int

        Static methods:
        global_serial_id_cnt    (property getter)
            global_serial_id_cnt -> int   (property getter)
        MEDPY_SET_global_serial_id_cnt(*args)
            MEDPY_SET_global_serial_id_cnt(int newval)

        Data descriptors:
            attributes -> StringFeatureAttrMapAdaptor   (property getter)
            data -> StringVecFloatMapAdaptor   (property getter)
        global_serial_id_cnt
            pid_pos_len -> IntPairIntIntMapAdaptor   (property getter)
            samples -> SampleVectorAdaptor   (property getter)
            tags -> StringUOSetStringMapAdaptor   (property getter)
            time_unit -> int   (property getter)
            weights    (property getter)


    class IdSamples
        Methods:
            id -> int   (property getter)
            samples -> SampleVectorAdaptor   (property getter)
            split -> int   (property getter)
            id <- int   (property setter)
            split <- int   (property setter)
            __init__(int _id) -> IdSamples
            __init__() -> IdSamples
            same_as(IdSamples other, int mode) -> bool
            set_split(int _split)

        Data descriptors:
            id -> int   (property getter)
            samples -> SampleVectorAdaptor   (property getter)
            split -> int   (property getter)


    class Mat
        Methods:
            avg    (property getter)
            missing_value -> float   (property getter)
            ncols -> int   (property getter)
            normalized_flag -> int   (property getter)
            nrows -> int   (property getter)
            row_ids    (property getter)
            signals -> StringVector   (property getter)
            std    (property getter)
            transposed_flag -> int   (property getter)
            missing_value <- float   (property setter)
            normalized_flag <- int   (property setter)
            row_ids <- int   (property setter)
            transposed_flag <- int   (property setter)
            __getitem__(IntVector index) -> float
            __init__() -> Mat
            __init__(int n_rows, int n_cols) -> Mat
            __init__(float * IN_ARRAY2) -> Mat
            __len__() -> unsigned long long
            __setitem__(IntVector index, float val)
            add_cols(Mat m_add)
            add_cols(float * m_add)
            add_rows(Mat m_add)
            add_rows(float * m_add)
            clear()
            get_col(int i_col)
            get_cols_avg_std()
            get_numpy_copy()
            get_numpy_view_unsafe()
            get_row(int i_row)
            get_sub_mat(vector< int > & rows_to_take, vector< int > & cols_to_take)
            get_sub_mat_by_flags(vector< int > & rows_to_take_flag, vector< int > & cols_to_take_flag)
            is_valid(bool output=False) -> bool
            is_valid() -> bool
            load(float * IN_ARRAY2)
            load(Mat x)
            load_numpy(float * IN_ARRAY2)
            load_transposed(float * IN_ARRAY2)
            normalize(int norm_type, float * wgts)
            normalize(int norm_type=Normalize_Cols)
            normalize()
            normalize(float * external_avg, float * external_std, int norm_type=1)
            normalize(float * external_avg, float * external_std)
            read_from_bin_file(string const & fname) -> int
            read_from_csv_file(string const & fname, int titles_line_flag) -> int
            resize(int n_rows, int n_cols)
            set_signals(StringVector sigs)
            set_val(float val)
            test(int n_rows, int n_cols)
            transpose()
            write_to_bin_file(string const & fname) -> int
            write_to_csv_file(string const & fname) -> int
            zero()

        Static methods:
        Normalize_Cols    (property getter)
            Normalize_Cols -> int   (property getter)
        Normalize_Rows    (property getter)
            Normalize_Rows -> int   (property getter)

        Data descriptors:
        Normalize_Cols
        Normalize_Rows
            avg    (property getter)
            missing_value -> float   (property getter)
            ncols -> int   (property getter)
            normalized_flag -> int   (property getter)
            nrows -> int   (property getter)
            row_ids    (property getter)
            signals -> StringVector   (property getter)
            std    (property getter)
            transposed_flag -> int   (property getter)


    class Model
        Methods:
            features -> Features   (property getter)
            verbosity -> int   (property getter)
            verbosity <- int   (property setter)
            __init__() -> Model
            add_age()
            add_feature_generator(string & name, string & signal)
            add_feature_generator_to_set(int i_set, string const & init_string)
            add_feature_generators(string & name, vector< string > & signals)
            add_feature_generators(string & name, vector< string > & signals, string init_string)
            add_feature_generators(string & name, string & signal, string init_string)
            add_feature_processor_to_set(int i_set, int duplicate, string const & init_string)
            add_gender()
            add_imputers()
            add_imputers(string init_string)
            add_imputers(vector< string > & features)
            add_imputers(vector< string > & features, string init_string)
            add_normalizers()
            add_normalizers(string init_string)
            add_normalizers(vector< string > & features)
            add_normalizers(vector< string > & features, string init_string)
            add_process_to_set(int i_set, int duplicate, string const & init_string)
            add_process_to_set(int i_set, string const & init_string)
            add_rep_processor_to_set(int i_set, string const & init_string)
            apply(PidRepository rep, Samples samples) -> int
            apply(PidRepository rep, Samples samples, int start_stage, int end_stage) -> int
            apply_feature_processors(Features features) -> int
            clear()
            collect_and_add_virtual_signals(PidRepository rep) -> int
            dprint_process(string const & pref, int rp_flag, int fg_flag, int fp_flag)
            filter_rep_processors()
            generate_all_features(PidRepository rep, Samples samples, Features features, StringVector req_feature_generators) -> int
            get_all_features_names(vector< string > & feat_names, int before_process_set)
            get_required_signal_names() -> StringVector
            init_from_json_file(std::string const & fname)
            init_from_json_file_with_alterations(std::string const & fname, StringVector json_alt) -> StringVector
            learn(PidRepository rep, Samples samples) -> int
            learn(PidRepository rep, Samples samples, int start_stage, int end_stage) -> int
            learn_and_apply_feature_processors(Features features) -> int
            learn_feature_generators(PidRepository rep, Samples learn_samples) -> int
            learn_feature_processors(Features features) -> int
            learn_rep_processors(PidRepository rep, Samples samples) -> int
            quick_learn_rep_processors(PidRepository rep, Samples samples) -> int
            read_from_file(string const & fname) -> int
            set_predictor(string name)
            set_predictor(string name, string init_string)
            write_feature_matrix(string const mat_fname) -> int
            write_to_file(std::string const & fname) -> int

        Data descriptors:
            features -> Features   (property getter)
            verbosity -> int   (property getter)


    class ModelStage
        Methods:
            __init__() -> ModelStage
        Data and other attributes:
        APPLY_FTR_GENERATORS = 2
        APPLY_FTR_PROCESSORS = 4
        APPLY_PREDICTOR = 6
        END = 8
        INSERT_PREDS = 7
        LEARN_FTR_GENERATORS = 1
        LEARN_FTR_PROCESSORS = 3
        LEARN_PREDICTOR = 5
        LEARN_REP_PROCESSORS = 0



    class PidRepository
        Methods:
            pids -> IntVector   (property getter)
            __init__() -> PidRepository
            dict_name(int_section_id, int_id) -> string
              returns name of section + id
            dict_prep_sets_lookup_table(int_section_id, list_String set_names) -> BoolVector
              returns a look-up-table for given set names
            dict_section_id(str_secName) -> int
              returns section id number for a given section name
            export_to_numpy(str_signame) -> SigExporter
              Returns the signal data represented as a list of numpy arrays, one for each field
        get_sig = __export_to_pandas(self, sig_name_str, translate=True, pids=None)
            get_sig(signame [, translate=True][, pids=None]) -> Pandas DataFrame
            translate : If True, will decode categorical fields into a readable representation in Pandas
            pid : If list is provided, will load only pids from the given list
                  If 'All' is provided, will use all available pids
            init(conf_file_name) -> int
            returns -1 if fails
            loadsig(str_signame) -> int
              load a signal
            read_all(conf_file_fname, [pids_to_take_array], [list_str_signals_to_take]) -> int
            returns -1 if fails
            reading a repository for a group of pids and signals.Empty group means all of it.
            read_all(conf_file_fname, [pids_to_take_array], [list_str_signals_to_take]) -> int
            returns -1 if fails
            reading a repository for a group of pids and signals.Empty group means all of it.
            read_all_i(std::string const & conf_fname, IntVector pids_to_take, IntVector signals_to_take) -> int
            sig_id(str_signame) -> int
              returns signal id number for a given signal name
            sig_type(str_signame) -> int
              returns signal type id for a given signal name
            uget(int_pid, int_sid) -> SigVectorAdaptor
              returns a vector of universal signals

        Data descriptors:
        dict
            PidRepository_dict_get() -> Dictionary
            pids -> IntVector   (property getter)


    class Sample
        Methods:
            id -> int   (property getter)
            outcome -> int   (property getter)
            outcomeTime -> int   (property getter)
            prediction    (property getter)
            split -> int   (property getter)
            time -> int   (property getter)
            id <- int   (property setter)
            outcome <- int   (property setter)
            outcomeTime <- int   (property setter)
            prediction <- float   (property setter)
            split <- int   (property setter)
            time <- int   (property setter)
            __copy__() -> Sample
            __init__() -> Sample
            parse_from_string(string & s, int time_unit) -> int
            print_(string const prefix)
            print_()
            write_to_string(string & s, int time_unit)

        Data descriptors:
            id -> int   (property getter)
            outcome -> int   (property getter)
            outcomeTime -> int   (property getter)
            prediction    (property getter)
            split -> int   (property getter)
            time -> int   (property getter)


    class Samples
        Methods:
            idSamples -> IdSamplesVectorAdaptor   (property getter)
            time_unit -> int   (property getter)
            time_unit <- int   (property setter)
            from_df(PandasAdaptor pandas_df)
            to_df() -> PandasAdaptor
            __init__() -> Samples
            append(Samples newSamples)
        as_df = __sample_export_to_pandas(self)
            clear()
            dilute(float prob)
            export_to_pandas_df() -> SampleVecExporter
            export_to_sample_vec() -> SampleVectorAdaptor
            get_attributes() -> StringVector
            get_categs()
            get_ids()
            get_predictions_size() -> int
            get_preds()
            get_str_attributes() -> StringVector
            get_y()
            insertRec(int pid, int time, float outcome, int outcomeTime)
            insertRec(int pid, int time, float outcome, int outcomeTime, float pred)
            insertRec(int pid, int time)
            insert_preds(Features featuresData) -> int
            nSamples() -> int
            normalize()
            read_from_bin_file(string const & file_name) -> int
            read_from_file(string const & file_name) -> int
            same_as(Samples other, int mode) -> bool
            sort_by_id_date()
            version() -> int
            write_to_bin_file(string const & file_name) -> int
            write_to_file(string const & fname) -> int

        Data descriptors:
            idSamples -> IdSamplesVectorAdaptor   (property getter)
            time_unit -> int   (property getter)


    class Sig
        Methods:
            __init__(Sig other) -> Sig
            date(int chan=0) -> int
            date() -> int
            days(int chan=0) -> int
            days() -> int
            hours(int chan=0) -> int
            hours() -> int
            minutes(int chan=0) -> int
            minutes() -> int
            months(int chan=0) -> int
            months() -> int
            time(int chan=0) -> int
            time() -> int
            timeU(int to_time_unit) -> int
            val(int chan=0) -> float
            val() -> float
            years(int chan=0) -> int
            years() -> int


    class Split
        Methods:
            nsplits -> int   (property getter)
            pid2split -> IntIntMapAdaptor   (property getter)
            __init__() -> Split
            clear()
            read_from_file(string const & fname) -> int
            write_to_file(string const & fname) -> int

        Data descriptors:
            nsplits -> int   (property getter)
            pid2split -> IntIntMapAdaptor   (property getter)


    class Time
        Methods:
            __init__() -> Time
        Data and other attributes:
        Date = 1
        DateTimeString = 7
        Days = 4
        Hours = 5
        Minutes = 6
        Months = 3
        Undefined = 0
        Years = 2

"
%enddef
