import { useState } from 'react'

export default function Tooltip({ term, children }) {
  const [visible, setVisible] = useState(false)

  return (
    <span className="tooltip-wrap">
      <span
        className="tooltip-trigger"
        onMouseEnter={() => setVisible(true)}
        onMouseLeave={() => setVisible(false)}
        onFocus={() => setVisible(true)}
        onBlur={() => setVisible(false)}
        tabIndex={0}
        role="button"
        aria-describedby="tooltip-desc"
      >
        {children}
      </span>
      {visible && (
        <span className="tooltip-box" role="tooltip" id="tooltip-desc">
          {term}
        </span>
      )}
    </span>
  )
}
