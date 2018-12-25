%define MPDOCSTRING
"


    class Dictionary
        Methods:
            __init__(MPDictionary self, PidRepository rep) -> Dictionary
            get_members_to_all_sets(Dictionary self, int section_id, int * members_array) -> IntVecIntMapAdaptor
            id(Dictionary self, string & signame) -> int
            name(Dictionary self, int id) -> string
            prep_sets_lookup_table(Dictionary self, int section_id, StringVector set_names)
            section_id(Dictionary self, string name) -> int

        Data descriptors:
        o

        Data and other attributes:


    class FeatureAttr
        Methods:
            imputed -> bool   (property getter)
            normalized -> bool   (property getter)
            imputed <- bool   (property setter)
            normalized <- bool   (property setter)
            __init__(MPFeatureAttr self) -> FeatureAttr

        Data descriptors:
        imputed
            imputed -> bool   (property getter)
        normalized
            normalized -> bool   (property getter)

        Data and other attributes:


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
            __init__(MPFeatures self) -> Features
            __init__(MPFeatures self, int _time_unit) -> Features
            __init__(MPFeatures self, Features other) -> Features
            append_samples(Features self, IdSamples in_samples)
            append_samples(Features self, Samples in_samples)
            clear(Features self)
            filter(Features self, StringVector selectedFeatures) -> int
        from_df = __features__from_df_imp(self, features_df)
            get_as_matrix(Features self, Mat mat)
            get_as_matrix(Features self, Mat mat, vector< string > names)
            get_as_matrix(Features self, Mat mat, vector< string > const names, int * int_in_buf)
            get_crc(Features self) -> unsigned int
            get_feature_names(Features self) -> StringVector
            get_max_serial_id_cnt(Features self) -> int
            get_pid_len(Features self, int pid) -> int
            get_pid_pos(Features self, int pid) -> int
            get_samples(Features self, Samples outSamples)
            init_all_samples(Features self, IdSamplesVectorAdaptor in_samples)
            init_pid_pos_len(Features self)
            insert_samples(Features self, IdSamples in_samples, int index)
            print_csv(Features self)
            read_from_csv_mat(Features self, string const & csv_fname) -> int
            set_as_matrix(Features self, Mat mat)
            set_time_unit(Features self, int _time_unit)
        to_df = __features__to_df_imp(self)
            version(Features self) -> int
            write_as_csv_mat(Features self, string const & csv_fname) -> int

        Static methods:
        global_serial_id_cnt    (property getter)
            global_serial_id_cnt -> int   (property getter)
        MEDPY_SET_global_serial_id_cnt(*args)
            MEDPY_SET_global_serial_id_cnt(int newval)

        Data descriptors:
        attributes
            attributes -> StringFeatureAttrMapAdaptor   (property getter)
        data
            data -> StringVecFloatMapAdaptor   (property getter)
        global_serial_id_cnt
        pid_pos_len
            pid_pos_len -> IntPairIntIntMapAdaptor   (property getter)
        samples
            samples -> SampleVectorAdaptor   (property getter)
        tags
            tags -> StringUOSetStringMapAdaptor   (property getter)
        time_unit
            time_unit -> int   (property getter)
        weights
            weights    (property getter)

        Data and other attributes:


    class IdSamples
        Methods:
            id -> int   (property getter)
            samples -> SampleVectorAdaptor   (property getter)
            split -> int   (property getter)
            id <- int   (property setter)
            split <- int   (property setter)
            __init__(MPIdSamples self, int _id) -> IdSamples
            __init__(MPIdSamples self) -> IdSamples
            same_as(IdSamples self, IdSamples other, int mode) -> bool
            set_split(IdSamples self, int _split)

        Data descriptors:
        id
            id -> int   (property getter)
        samples
            samples -> SampleVectorAdaptor   (property getter)
        split
            split -> int   (property getter)

        Data and other attributes:


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
            __getitem__(Mat self, IntVector index) -> float
            __init__(MPMat self) -> Mat
            __init__(MPMat self, int n_rows, int n_cols) -> Mat
            __init__(MPMat self, float * IN_ARRAY2) -> Mat
            __len__(Mat self) -> unsigned long long
            __setitem__(Mat self, IntVector index, float val)
            add_cols(Mat self, Mat m_add)
            add_cols(Mat self, float * m_add)
            add_rows(Mat self, Mat m_add)
            add_rows(Mat self, float * m_add)
            clear(Mat self)
            get_col(Mat self, int i_col)
            get_cols_avg_std(Mat self)
            get_numpy_copy(Mat self)
            get_numpy_view_unsafe(Mat self)
            get_row(Mat self, int i_row)
            get_sub_mat(Mat self, vector< int > & rows_to_take, vector< int > & cols_to_take)
            get_sub_mat_by_flags(Mat self, vector< int > & rows_to_take_flag, vector< int > & cols_to_take_flag)
            is_valid(Mat self, bool output=False) -> bool
            is_valid(Mat self) -> bool
            load(Mat self, float * IN_ARRAY2)
            load(Mat self, Mat x)
            load_numpy(Mat self, float * IN_ARRAY2)
            load_transposed(Mat self, float * IN_ARRAY2)
            normalize(Mat self, int norm_type, float * wgts)
            normalize(Mat self, int norm_type=Normalize_Cols)
            normalize(Mat self)
            normalize(Mat self, float * external_avg, float * external_std, int norm_type=1)
            normalize(Mat self, float * external_avg, float * external_std)
            read_from_bin_file(Mat self, string const & fname) -> int
            read_from_csv_file(Mat self, string const & fname, int titles_line_flag) -> int
            resize(Mat self, int n_rows, int n_cols)
            set_signals(Mat self, StringVector sigs)
            set_val(Mat self, float val)
            test(Mat self, int n_rows, int n_cols)
            transpose(Mat self)
            write_to_bin_file(Mat self, string const & fname) -> int
            write_to_csv_file(Mat self, string const & fname) -> int
            zero(Mat self)

        Static methods:
        Normalize_Cols    (property getter)
            Normalize_Cols -> int   (property getter)
        Normalize_Rows    (property getter)
            Normalize_Rows -> int   (property getter)

        Data descriptors:
        Normalize_Cols
        Normalize_Rows
        avg
            avg    (property getter)
        missing_value
            missing_value -> float   (property getter)
        ncols
            ncols -> int   (property getter)
        normalized_flag
            normalized_flag -> int   (property getter)
        nrows
            nrows -> int   (property getter)
        row_ids
            row_ids    (property getter)
        signals
            signals -> StringVector   (property getter)
        std
            std    (property getter)
        transposed_flag
            transposed_flag -> int   (property getter)

        Data and other attributes:


    class Model
        Methods:
            features -> Features   (property getter)
            verbosity -> int   (property getter)
            verbosity <- int   (property setter)
            __init__(MPModel self) -> Model
            add_age(Model self)
            add_feature_generator(Model self, string & name, string & signal)
            add_feature_generator_to_set(Model self, int i_set, string const & init_string)
            add_feature_generators(Model self, string & name, vector< string > & signals)
            add_feature_generators(Model self, string & name, vector< string > & signals, string init_string)
            add_feature_generators(Model self, string & name, string & signal, string init_string)
            add_feature_processor_to_set(Model self, int i_set, int duplicate, string const & init_string)
            add_gender(Model self)
            add_imputers(Model self)
            add_imputers(Model self, string init_string)
            add_imputers(Model self, vector< string > & features)
            add_imputers(Model self, vector< string > & features, string init_string)
            add_normalizers(Model self)
            add_normalizers(Model self, string init_string)
            add_normalizers(Model self, vector< string > & features)
            add_normalizers(Model self, vector< string > & features, string init_string)
            add_process_to_set(Model self, int i_set, int duplicate, string const & init_string)
            add_process_to_set(Model self, int i_set, string const & init_string)
            add_rep_processor_to_set(Model self, int i_set, string const & init_string)
            apply(Model self, PidRepository rep, Samples samples) -> int
            apply(Model self, PidRepository rep, Samples samples, int start_stage, int end_stage) -> int
            apply_feature_processors(Model self, Features features) -> int
            clear(Model self)
            collect_and_add_virtual_signals(Model self, PidRepository rep) -> int
            dprint_process(Model self, string const & pref, int rp_flag, int fg_flag, int fp_flag)
            filter_rep_processors(Model self)
            generate_all_features(Model self, PidRepository rep, Samples samples, Features features, StringVector req_feature_generators) -> int
            get_all_features_names(Model self, vector< string > & feat_names, int before_process_set)
            get_required_signal_names(Model self) -> StringVector
            init_from_json_file(Model self, std::string const & fname)
            init_from_json_file_with_alterations(Model self, std::string const & fname, StringVector json_alt) -> StringVector
            learn(Model self, PidRepository rep, Samples samples) -> int
            learn(Model self, PidRepository rep, Samples samples, int start_stage, int end_stage) -> int
            learn_and_apply_feature_processors(Model self, Features features) -> int
            learn_feature_generators(Model self, PidRepository rep, Samples learn_samples) -> int
            learn_feature_processors(Model self, Features features) -> int
            learn_rep_processors(Model self, PidRepository rep, Samples samples) -> int
            quick_learn_rep_processors(Model self, PidRepository rep, Samples samples) -> int
            read_from_file(Model self, string const & fname) -> int
            set_predictor(Model self, string name)
            set_predictor(Model self, string name, string init_string)
            write_feature_matrix(Model self, string const mat_fname) -> int
            write_to_file(Model self, std::string const & fname) -> int

        Data descriptors:
        features
            features -> Features   (property getter)
        verbosity
            verbosity -> int   (property getter)

        Data and other attributes:


    class ModelStage
        Methods:
            __init__(MPModelStage self) -> ModelStage

        Data descriptors:

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
            __init__(MPPidRepository self) -> PidRepository
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
            read_all_i(PidRepository self, std::string const & conf_fname, IntVector pids_to_take, IntVector signals_to_take) -> int
            sig_id(str_signame) -> int
              returns signal id number for a given signal name
            sig_type(str_signame) -> int
              returns signal type id for a given signal name
            uget(int_pid, int_sid) -> SigVectorAdaptor
              returns a vector of universal signals

        Data descriptors:
        dict
            PidRepository_dict_get(PidRepository self) -> Dictionary
        o
        pids
            pids -> IntVector   (property getter)

        Data and other attributes:


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
            __copy__(Sample self) -> Sample
            __init__(MPSample self) -> Sample
            parse_from_string(Sample self, string & s, int time_unit) -> int
            print_(Sample self, string const prefix)
            print_(Sample self)
            write_to_string(Sample self, string & s, int time_unit)

        Data descriptors:
        id
            id -> int   (property getter)
        outcome
            outcome -> int   (property getter)
        outcomeTime
            outcomeTime -> int   (property getter)
        prediction
            prediction    (property getter)
        split
            split -> int   (property getter)
        time
            time -> int   (property getter)

        Data and other attributes:


    class Samples
        Methods:
            idSamples -> IdSamplesVectorAdaptor   (property getter)
            time_unit -> int   (property getter)
            time_unit <- int   (property setter)
            MEDPY__from_df(Samples self, PandasAdaptor pandas_df)
            MEDPY__to_df(Samples self) -> PandasAdaptor
            __init__(MPSamples self) -> Samples
            append(Samples self, Samples newSamples)
        as_df = __sample_export_to_pandas(self)
            clear(Samples self)
            dilute(Samples self, float prob)
            export_to_pandas_df(Samples self) -> SampleVecExporter
            export_to_sample_vec(Samples self) -> SampleVectorAdaptor
        from_df = __from_df_imp(self, df)
            get_attributes(Samples self) -> StringVector
            get_categs(Samples self)
            get_ids(Samples self)
            get_predictions_size(Samples self) -> int
            get_preds(Samples self)
            get_str_attributes(Samples self) -> StringVector
            get_y(Samples self)
            insertRec(Samples self, int pid, int time, float outcome, int outcomeTime)
            insertRec(Samples self, int pid, int time, float outcome, int outcomeTime, float pred)
            insertRec(Samples self, int pid, int time)
            insert_preds(Samples self, Features featuresData) -> int
            nSamples(Samples self) -> int
            normalize(Samples self)
            read_from_bin_file(Samples self, string const & file_name) -> int
            read_from_file(Samples self, string const & file_name) -> int
            same_as(Samples self, Samples other, int mode) -> bool
            sort_by_id_date(Samples self)
        to_df = __to_df_imp(self)
            version(Samples self) -> int
            write_to_bin_file(Samples self, string const & file_name) -> int
            write_to_file(Samples self, string const & fname) -> int

        Data descriptors:
        idSamples
            idSamples -> IdSamplesVectorAdaptor   (property getter)
        time_unit
            time_unit -> int   (property getter)

        Data and other attributes:


    class Sig
        Methods:
            __init__(MPSig self, Sig other) -> Sig
            date(Sig self, int chan=0) -> int
            date(Sig self) -> int
            days(Sig self, int chan=0) -> int
            days(Sig self) -> int
            hours(Sig self, int chan=0) -> int
            hours(Sig self) -> int
            minutes(Sig self, int chan=0) -> int
            minutes(Sig self) -> int
            months(Sig self, int chan=0) -> int
            months(Sig self) -> int
            time(Sig self, int chan=0) -> int
            time(Sig self) -> int
            timeU(Sig self, int to_time_unit) -> int
            val(Sig self, int chan=0) -> float
            val(Sig self) -> float
            years(Sig self, int chan=0) -> int
            years(Sig self) -> int

        Data descriptors:

        Data and other attributes:


    class Split
        Methods:
            nsplits -> int   (property getter)
            pid2split -> IntIntMapAdaptor   (property getter)
            __init__(MPSplit self) -> Split
            clear(Split self)
            read_from_file(Split self, string const & fname) -> int
            write_to_file(Split self, string const & fname) -> int

        Data descriptors:
        nsplits
            nsplits -> int   (property getter)
        o
        pid2split
            pid2split -> IntIntMapAdaptor   (property getter)

        Data and other attributes:


    class Time
        Methods:
            __init__(MPTime self) -> Time

        Data descriptors:

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
