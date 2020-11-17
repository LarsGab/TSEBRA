class Features:
    def __init__(self):
        self.name2class = {'numb_introns' : Float_feature, \
            'transcript_length' : Float_feature, 'intron_length' : Float_feature, \
            'fraction_intron_leng' : Float_feature}
        self.features = {}

    def add(self, f_name, value, pref = None):
        arg = [value]
        if not pref == None:
            arg.append(pref)
        if f_name in self.features.keys():
            self.features[f_name] = self.name2class[f_name](*arg)
        elif f_name in self.feature_names():
            self.features.update({f_name : self.name2class[f_name](*arg)})
        else:
            eprint("Can't add {}. Feature not available.".format(f_name))

    def get_value(self, f_name):
        if f_name in self.features.keys():
            return self.features[f_name].value
        elif f_name in self.feature_names():
            eprint("Can't get {}. Feature has not been added yet.".format(f_name))
        else:
            eprint("Can't get {}. Feature not available.".format(f_name))

    def get_value_list(self, f_names = []):
        if not f_names:
            f_names = self.features.keys()

        result = []
        for k in f_names:
            result.append(self.features[k].value)
        return result

    def get_dist_list(self, f_names = []):
        if not f_names:
            f_names = self.features.keys()
        result = []
        for k in f_names:
            result.append(self.features[k].dist())
        return result

    def feature_names(self):
        return self.name2class.keys()


class Float_feature:
    def __init__(self, value, pref=0):
        self.value = value
        self.pref = pref

    def dist(self):
        return abs(self.value - self.pref)

class Intron_numb:
    def __init__(self, i_numb, pref = 2):
        self.intron_numb = i_numb
        self.pref = pref

    def dist(self):
        return abs(self.intron_numb - pref)

class Transcript_length:
    def __init__(self, length, pref = 0):
        self.length = length
        self.pref = pref

    def dist(self):
        return abs(self.length - pref)

class Intron_length:
    def __init__(self, length, pref = 0):
        self.length = length
        self.pref = pref

    def dist(self):
        return abs(self.length - pref)


class Intron_length_fraction:
        def __init__(self, fraction, pref = 1):
            self.fraction = fraction
            self.pref = pref

        def dist(self):
            return abs(self.fraction - pref)
