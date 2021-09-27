### CLASS ----------------------------------------------------------------------
student <- setClass("student", 
                    slots=list(name="character", 
                               age="numeric", 
                               GPA="numeric"))
### METHODS --------------------------------------------------------------------
setMethod("show",
          "student",
          function(object) {
            cat(object@name, "\n")
            cat(object@age, "years old\n")
            cat("GPA:", object@GPA, "\n")})

setMethod("increaseAge",
          "student",
          function(object) {
            cat(object@name, "\n")
            cat(object@age, "years old\n")
            cat("GPA:", object@GPA, "\n")})

### OBJECTS --------------------------------------------------------------------
student_1 <- student(name="John", age=21, GPA=3.5)