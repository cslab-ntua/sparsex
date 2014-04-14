/*
 * Encodings.hpp -- All about encoding types
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef ENCODINGS_HPP
#define ENCODINGS_HPP

#include <boost/bimap.hpp>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

namespace csx {

class Encoding
{
public:
    enum Type {
        None,
        Horizontal,
        Vertical,
        Diagonal,
        AntiDiagonal,
        BlockRow1,
        BlockRow2,
        BlockRow3,
        BlockRow4,
        BlockRow5,
        BlockRow6,
        BlockRow7,
        BlockRow8,
        BlockCol1,
        BlockCol2,
        BlockCol3,
        BlockCol4,
        BlockCol5,
        BlockCol6,
        BlockCol7,
        BlockCol8,
        __EndOfTypes__,
        // Group definitions
        BlockRows,
        BlockCols,
        All,
        __EndOfGroups__,
    };

    typedef pair<Type, size_t> Instantiation;

    Encoding(Type type)
        : type_(type)
    { }

    Encoding(const string &name)
        : type_(rev_names_[name])
    { }

    Type GetType() const
    {
        return type_;
    }

    const string &GetFullName() const
    {
        return names_[type_].second;
    }

    const string &GetShortName() const
    {
        return names_[type_].first;
    }

    bool IsGroup() const
    {
        return type_ > __EndOfTypes__;
    }

    bool IsBlockRow() const
    {
        return type_ >= BlockRowMin && type_ <= BlockRowMax;
    }
    
    bool IsBlockCol() const
    {
        return type_ >= BlockColMin && type_ <= BlockColMax;
    }

    bool IsBlock() const
    {
        return IsBlockRow() || IsBlockCol();
    }


    size_t GetBlockAlignment() const
    {
        size_t ret;
        if (IsBlockRow())
            ret = type_ - BlockRowMin + 1;
        else if (IsBlockCol())
            ret = type_ - BlockColMin + 1;
        else
            ret = 0;
        return ret;
    }

    const void GetTypes(vector<Type> &types) const;

    static const Type Max = __EndOfTypes__;
    static const Type BlockRowMin = BlockRow1;
    static const Type BlockRowMax = BlockRow8;
    static const Type BlockColMin = BlockCol1;
    static const Type BlockColMax = BlockCol8;
private:
    static std::map<string, Type> InitInverseNameMap();

    Type type_;
    static std::map<Type, pair<string, string> > names_;
    static std::map<string, Type> rev_names_;
};

ostream &operator<<(ostream &out, Encoding &enc)
{
    out << enc.GetFullName();
    return out;
}

Encoding::Type &operator++(Encoding::Type &type)
{
    assert(type < Encoding::Max);
    type = static_cast<Encoding::Type>(type+1);
    return type;
}

Encoding::Type operator++(Encoding::Type &type, int)
{
    assert(type < Encoding::Max);
    Encoding::Type orig = type;
    type = static_cast<Encoding::Type>(type+1);
    return orig;
}

class EncodingSequence
{
public:
    typedef vector<pair<Encoding, vector<size_t> > >::const_iterator
    const_iterator;
    typedef vector<size_t>::const_iterator delta_const_iterator;

    EncodingSequence(const string &seqstr);

    size_t GetSize() const
    {
        return sequence_.size();
    }

    const_iterator Cbegin() const
    {
        return sequence_.begin();
    }
    
    const_iterator Cend() const
    {
        return sequence_.end();
    }

    delta_const_iterator DeltasCbegin(const_iterator &iter) const
    {
        const vector<size_t> &deltas = iter->second;
        return deltas.begin();
    }

    delta_const_iterator DeltasCend(const_iterator &iter) const
    {
        const vector<size_t> &deltas = iter->second;
        return deltas.end();
    }

    bool IsExplicit() const
    {
        return explicit_;
    }

    void Print(ostream &out) const;

private:
    vector<pair<Encoding, vector<size_t> > > sequence_;
    bool explicit_;
};

ostream &operator<<(ostream &out, EncodingSequence &encseq)
{
    encseq.Print(out);
    return out;
}

class PreprocessingMethod
{
public:
    enum Type {
        None,
        FixedWindow,
        FixedPortion,
    };

    PreprocessingMethod(Type meth)
        : method_type_(meth)
    {
        InitMethodNames();
    }

    PreprocessingMethod(const string &name)
    {
        InitMethodNames();
        method_type_ = method_names_.right.at(name);
    }

    const string &GetName() const
    {
        return method_names_.left.at(method_type_);
    }

    Type GetType() const
    {
        return method_type_;
    }

private:
    static void InitMethodNames()
    {
        if (method_names_.empty()) {
            method_names_.insert(NameMap::value_type(None, "none"));
            method_names_.insert(NameMap::value_type(FixedWindow, "window"));
            method_names_.insert(NameMap::value_type(FixedPortion, "portion"));
        }
    }

    typedef boost::bimap<Type, string> NameMap;
    static NameMap method_names_;
    Type method_type_;
};

ostream &operator<<(ostream &out, PreprocessingMethod &meth)
{
    out << meth.GetName();
    return out;
}

}   // end of namespace csx

#endif  // ENCODINGS_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
