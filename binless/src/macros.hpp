#ifndef MACROS_HPP
#define MACROS_HPP

#define SQUARE(x) ((x)*(x))

//define a private variable and a public getter / setter pair
//CType is the type of the variable as passed in the setter (e.g. const Type&)
#define BINLESS_GET_SET_DECL(Type, CType, Name)                  \
public:                                                          \
  Type get_##Name() const {                                      \
    return Name##_;                                              \
  };                                                             \
  void set_##Name(CType value) {                                 \
    Name##_ = value;                                             \
  }                                                              \
  private:                                                       \
    Type Name##_;
  
  //define a const private variable and public getter
#define BINLESS_GET_CONST_DECL(Type, Name)                         \
  public:                                                          \
    Type get_##Name() const {                                      \
      return Name##_;                                              \
    };                                                             \
    private:                                                       \
      const Type Name##_;
    
    //define a private variable and public getter that returns a reference to it
#define BINLESS_GET_REF_DECL(Type, Name)                             \
    public:                                                          \
      Type& get_##Name() { return Name##_; }                         \
      const Type& get_##Name() const { return Name##_; }             \
      private:                                                       \
        Type Name##_;

#endif

