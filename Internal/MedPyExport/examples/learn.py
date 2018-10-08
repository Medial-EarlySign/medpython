
# coding: utf-8

# In[1]:

import sys
import numpy as np
import pandas as pd
#sys.path.insert(0,'/nas1/UsersData/shlomi/MR/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release/MedPython')

import medpython as med


# In[2]:

class conf:
    NAME='config for dneph full model'
    #STEPS  Train
    STEPS=['Train','Validate']
    #STEPS  Validate
    
    REPOSITORY='/home/Repositories/THIN/thin_jun2017/thin.repository'

    SPLITS='/home/Repositories/THIN/thin_jun2017/splits6.thin'

    # cmd line override with --train_test_learning_samples_fname
    LEARNING_SAMPLES='/nas1/UsersData/shlomi/MR/for_shlomi/mi_train.samples'

    LEARNING_SAMPLES_PROB=0.25

    # cmd line override with --train_test_validation_samples_fname
    VALIDATION_SAMPLES='/nas1/UsersData/shlomi/MR/for_shlomi/mi_test.samples'

    # cmd line override with --work_dir
    WORK_DIR='/nas1/UsersData/shlomi/MR/for_shlomi/res_mi'
    #WORK_DIR='/home/python/for_shlomi/res_mi'

    # cmd line override with --train_test_prefix
    PREFIX='mi'

    VALIDATION_PREFIX='mi_test'

    # cmd line override with --train_test_model_init_fname
    MODEL_INIT_FILE='/nas1/UsersData/shlomi/MR/for_shlomi/mi_model.json'

    ACTIVE_SPLITS='full'
    #ACTIVE_SPLITS  0

    CV_ON_TRAIN1=0
    SAVE_MATRIX=0


# In[3]:

class TrainTest:
    
    def __init__(self, conf):
        self.conf = conf
        self.can_write_files = False
    
    def print_performance(self,s, prefix, split, strm):
        print 'TT: print_performance() not implemented'
    
    def print_distrib(self,s):
        y = s.get_y();
        freq = {}
        for i in range(len(y)):
            freq[y[i]] = freq.get(y[i],0) +1
        top_5 = sorted(list(freq.iteritems()),key=lambda i:i[1], reverse=True)[:5]
        print "TT: patients size: {} y size: {} distribution: {}".format(len(s.idSamples), len(y), top_5.__repr__())
        
    def read_all_input_files(self):
        m = med.Model()
        self.json_alt = []
        self.json_alt = m.init_from_json_file_with_alterations(self.conf.MODEL_INIT_FILE, self.json_alt)
        m.init_from_json_file(self.conf.MODEL_INIT_FILE)

        self.learning_samples = med.Samples()
        self.learning_samples.read_from_file(self.conf.LEARNING_SAMPLES)
        self.learning_samples.dilute(self.conf.LEARNING_SAMPLES_PROB)
        self.learning_samples.time_unit = med.Time.Date

        self.validation_samples = med.Samples()
        if self.conf.VALIDATION_SAMPLES:
            self.validation_samples.read_from_file(self.conf.VALIDATION_SAMPLES)
        self.validation_samples.time_unit = med.Time.Date

        self.sp = med.Split()
        self.sp.read_from_file(self.conf.SPLITS)

        self.active_splits={}
        for split in range(self.sp.nsplits):
            self.active_splits[split] = False         #do_cv_for_all_splits
        self.active_splits[self.sp.nsplits] = True    #do_full_model
        for split in range(self.sp.nsplits+1):
            print "TT: split [{}] active = {}".format(split, self.active_splits[split])

        req_names = m.get_required_signal_names()
        sigs = ["TRAIN"] + list(req_names)
        sigs = sorted(set(sigs))            # Sort+uniq

        pids = np.concatenate((self.learning_samples.get_ids(), 
                               self.validation_samples.get_ids()), axis=0)

        self.rep = med.PidRepository()
        if self.rep.read_all(self.conf.REPOSITORY, pids , sigs)<0:
            print "TT: ERROR could not read repository {}".format(self.conf.REPOSITORY);

        self.paramList = []

    def train(self):
        self.models = []
        v=[]
        l=[]
        for x in range(self.sp.nsplits+1):
            self.models.append(med.Model())
            v.append(med.Samples())
            l.append(med.Samples())


        learningCVPredictions = med.Samples()
        learningCVPredictions.clear()
        featuresVecLearn=[]
        featuresVecApply=[]


        for split in range(self.sp.nsplits+1):
            if split == self.sp.nsplits:
                print "TT: full model (denoted as split {})".format(split)
            else:
                print "TT: split {} (out of {})".format(split, self.sp.nsplits);
            if not self.active_splits[split]:
                print "TT: skipping split {} as it was not specified".format(split);
                continue
            med.Features.global_serial_id_cnt = 0
            print "learning_samples.idSamples.size() = {}".format(len(self.learning_samples.idSamples))

            for sample_i in range(len(self.learning_samples.idSamples)):
                sample = self.learning_samples.idSamples[sample_i]
                if len(self.sp.pid2split)>0: mysplit = self.sp.pid2split[sample.id]
                else: mysplit = sample.split
                if mysplit != split :
                    if len(l[split].idSamples) < 3:
                        print "TT: pid: {} belongs to split: {} = learning".format(sample.id, mysplit);
                    l[split].idSamples.append(sample);
                else:
                    if len(v[split].idSamples) < 3:
                        print "TT: pid: {} belongs to split: {} = validating".format(sample.id, mysplit);
                    v[split].idSamples.append(sample);
 
            self.models[split] = med.Model();

            my_json_alt = list(self.json_alt)
            my_json_alt.append('_SPLIT_::' + str(split))
            my_json_alt.append('_WORK_FOLDER_::' + self.conf.WORK_DIR)
 
            self.models[split].init_from_json_file_with_alterations(self.conf.MODEL_INIT_FILE, my_json_alt);

            print "TT: learning split {} on learning samples with distribution:".format(split)
            self.print_distrib(l[split]);

            if self.models[split].learn(self.rep, l[split])<0:
                print "TT: Learning model for split {} failed".format(split)

            print "TT: Collected CV predictions for [{}] patients on the learning samples".format(len(learningCVPredictions.idSamples))

            model_f = self.conf.WORK_DIR + "/" + self.conf.PREFIX + "_S" + str(split) + ".model"

            print "TT: split {} : writing model file {}".format(split, model_f);
            if self.models[split].write_to_file(model_f) < 0:
                print "TT: split {} : FAILED writing model file {}".format(split, model_f);
        
    def validate(self):
        print "TT: ********* validate *********"
        should_go_in_cv = set()
        #split_validate_to_cv_and_full(should_go_in_cv)
        if self.conf.CV_ON_TRAIN1==1:
            raise Exception("no implemented")
        else:
            #simply put all learning_samples in should_go_to_cv
            for orec_i in range(len(self.learning_samples.idSamples)):
                should_go_in_cv.add(self.learning_samples.idSamples[orec_i].id)
            print "TT: cv_on_train_flag=[{}], found [{}]/[{}] pids in validation_samples that will go into CV instead of full model".format(
                self.conf.CV_ON_TRAIN1, len(should_go_in_cv), len(self.validation_samples.idSamples))

        self.models = []
        for x in range(self.sp.nsplits+1):
            self.models.append(med.Model())

        validationCVPredictions = med.Samples()
        validationFullModelPredictions = med.Samples()

        if self.can_write_files:
            perf_file = self.conf.WORK_DIR + '/' + self.conf.VALIDATION_PREFIX + '_perf.csv'
            perf_fstream = open(perf_file, 'w')
            if not perf_fstream:
                print 'can''t open file {} for write'.format(perf_file)
            else: 
                print 'opening file {} for write'.format(perf_file)
                perf_fstream.write('prefix,split,auc,ctrls,cases,total,corr')
        else: perf_fstream = None

        for split in range(self.sp.nsplits+1):
            if split == self.sp.nsplits:
                print 'TT: Applying full model (denoted as split {}) on all validation samples that were not in the learning set'.format(
                split)
            else: print 'TT: Applying split {} model (out of {}) on validation samples that were in the learning set'.format(split, self.sp.nsplits)
            if not self.active_splits[split]:
                print 'TT: skipping {} as it was not specified'.format(split)
                continue
            v = med.Samples()

            for sample_i in range(len(self.validation_samples.idSamples)):
                sample = self.validation_samples.idSamples[sample_i]
                if split == self.sp.nsplits and not sample.id in should_go_in_cv:
                    v.idSamples.append(sample)
                if split < self.sp.nsplits and sample.id in should_go_in_cv and self.sp.pid2split[sample.id] == split:
                    v.idSamples.append(sample)
            if len(v.idSamples) == 0:
                print "TT: no samples were found for split {}, continuing...".format(split)
                continue

            print 'TT: split {} has [{}] patients '.format(split, len(v.idSamples))

            self.models[split] = med.Model()
            model_f = self.conf.WORK_DIR + '/' + self.conf.PREFIX + '_S' + str(split) + '.model'
            print 'TT: before reading model file {}'.format(model_f)
            if self.models[split].read_from_file(model_f) < 0:
                print 'TT: Validation: split {} : FAILED reading models file {}'.format(split, model_f)
                raise Exception('FAILED')
            print 'TT: completed reading model file {}'.format(model_f)
            if self.models[split].apply(self.rep, v) < 0:
                print 'TT: Applying model for split {} failed'.format(split)
                raise Exception('FAILED')
            print 'TT: After apply of split {} from file {}'.format(split, model_f)
            x = med.Mat()
            self.models[split].features.get_as_matrix(x)

            for sample_i in range(len(v.idSamples)):
                v.idSamples[sample_i].set_split(split)

            if split < self.sp.nsplits:
                print 'TT: split {} : predictions CV performance on validation samples that were used in train'.format(split)
                self.print_performance(v, validation_prefix, split, perf_fstream)
                validationCVPredictions.append(v);
            else:
                validationFullModelPredictions.append(v)

        if len(validationCVPredictions.idSamples) > 0:
            split = 999
            print 'TT: split {} : combined cross validated predictions performance on validation samples that were used in train'.format(split)
            self.print_performance(validationCVPredictions, self.conf.VALIDATION_PREFIX, split, perf_fstream)
        if len(validationFullModelPredictions.idSamples) > 0:
            split = self.sp.nsplits
            print 'TT: split {} : full model predictions performance on validation samples outside of train'.format(split)
            self.print_performance(validationFullModelPredictions, self.conf.VALIDATION_PREFIX, split, perf_fstream)

        validationCombinedPredictions = med.Samples()
        validationCombinedPredictions.append(validationCVPredictions)
        validationCombinedPredictions.append(validationFullModelPredictions)

        # combined performance
        if len(validationCombinedPredictions.idSamples) > 0:
            split = 999999
            print ('TT: split {} : combined cross validated + full model predictions performance on all samples'                    '(cv prediction for samples used in train, full model for samples not in train)'.format(split))
            self.print_performance(validationCombinedPredictions, self.conf.VALIDATION_PREFIX, split, perf_fstream);
        f_pl = self.conf.WORK_DIR+'/'+self.conf.VALIDATION_PREFIX+'.preds'
        validationCombinedPredictions.write_to_file(f_pl)
        if self.can_write_files:
            perf_fstream.close()
            print 'TT: wrote [{}]'.format(perf_file)
        


# In[4]:

tt = TrainTest(conf)
tt.read_all_input_files()
print med.cerr()


# In[5]:

tt.train()


# In[6]:

print med.cerr()


# In[7]:

tt.validate()


# In[8]:

print med.cerr()


# In[ ]:




