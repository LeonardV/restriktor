# Long output text is wrapped at getOption("restriktor.linewidth", 76) characters
# (see R/output_wrap.R). In the tests, this is switched off, such that the
# expected (error) messages can be matched as one line; the wrapping itself is
# tested in test-output-wrap-messages.R.
options(restriktor.linewidth = Inf)
