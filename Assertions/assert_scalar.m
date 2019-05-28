function assert_scalar(name, value)
  assert(numel(value) == 1, [name ' must be a scalar.']);
end

